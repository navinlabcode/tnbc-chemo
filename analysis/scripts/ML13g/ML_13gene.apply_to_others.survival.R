suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Apply the 13-gene model to other cohorts with survival data
#
# Simplifed from logistic_feature_selection.response_ARTEMIS.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025-03-12
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr(base_size = 8, legend = "right") %+replace% theme(axis.ticks.length = unit(0.1, "inch")))
library(cli)
library(tictoc)
library(glue)
library(scales)
library(tools)
library(survminer)
library(survival)
library(caret)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
#------ parse command line arguments ------
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    study_name <- cmdargs[1]
    f_type <- cmdargs[2]
} else {
    study_name <- "METABRIC"
    f_type <- "chemoYes"
    # study_name <- "SCANB"
    # f_type <- "dUTP"
}
#------ load the 13-gene model  ------

dir_model <- file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102",
    "lv01.aneuploidy_tri_type.aneuploid.pure5",
    "survival_analysis",
    "artemis_chemoyes.logistic_informative_genes",
    "cancer_mm_archetype_markers_CellstatePCRrelated_CellstateMarkers"
)

mulv_model <- read_rds(file.path(
    dir_model,
    sprintf("%s.multivariate_logistic_regression_model.rds", "pCR_status")
))
mulv_model_genes <- read_lines(file.path(
    dir_model,
    sprintf("%s.multivariate_logistic_regression_model_genes.txt", "pCR_status")
))
mulv_model_genes <- intersect(names(mulv_model$coefficients), mulv_model_genes)
#------ model gene matrix ------
mat_model <- read_rds(file.path("/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/psbulk/", "normalized.matrix.rds"))
#------ prepare the input expression data ------

# # Only survival data; no response data
if (study_name == "METABRIC") {
    # metabric's rest dataset
    f_mat <- file.path("/volumes/USR1/yyan/shared/data/METABRIC/ruli/", "tnbc_chemoYES_expression.mat.rds")
    f_df <- file.path("/volumes/USR1/yyan/shared/data/METABRIC/ruli/", "tnbc_chemoYES_survival.df.rds")
}

if (study_name == "SCANB") {
    # SCANB TNBC updated 2025-03-04
    if (f_type != "ALL_protocols") {
        dir_lib <- file.path("/volumes/USR1/yyan/shared/data/SCANB_TNBC_chemo/stringtie_fpkm_gene_data_unadjusted/", f_type)
    } else {
        dir_lib <- file.path("/volumes/USR1/yyan/shared/data/SCANB_TNBC_chemo/stringtie_fpkm_gene_data_LibProtocolAdjusted/", f_type)
    }
    f_df <- file.path(dir_lib, "scanb_tnbc.df_clinical.rds")
    f_mat <- file.path(dir_lib, "scanb_tnbc.expr_mat.rds")
}

dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas/summary_fig/ML13gene_app",
    study_name, "ml_13gene_survival", f_type
)
fs::dir_create(dir_res)

#------------------ ~~~ START ~~~ --------------------
mat <- read_rds(f_mat)
df <- read_rds(f_df)
mat <- as.matrix(mat)
head(df)
head(colnames(mat))

if (study_name == "METABRIC") {
    try(df$Patient.ID <- df$METABRIC_ID)
    df <- column_to_rownames(df, var = "Patient.ID")
    print(table(df$OS_STATUS))
    eval_x_choices <- c("OS_DAYS")
    eval_y_choices <- c("OS_STATUS")
}

if (study_name == "SCANB") {
    df <- as.data.frame(df)
    rownames(df) <- df$GEX.assay
    eval_x_choices <- c("DRFi_days", "OS_days", "RFi_days", "BCFi_days")
    eval_y_choices <- c("DRFi_event", "OS_event", "RFi_event", "BCFi_event")
    eval_x_choices <- c("OS_days")
    eval_y_choices <- c("OS_event")
}


if (study_name == "artemis_pretx") {
    rownames(df) <- df$patient
    df <- df[colnames(mat), ]
}

if (study_name == "ISPY990") {
    df <- as.data.frame(df)
    rownames(df) <- df[["Patient Identifier"]]
    mat[is.na(mat)] <- 0
}

shared_pat <- intersect(rownames(df), colnames(mat))
str(shared_pat)
df <- df[shared_pat, ]
mat <- mat[, shared_pat]

print(range(mat))

stopifnot(identical(rownames(df), colnames(mat)))

if (!all(mulv_model_genes %in% rownames(mat))) {
    cat("These genes are not present:\n")
    str(mulv_model_genes[!mulv_model_genes %in% rownames(mat)])
}


idx <- match(mulv_model_genes, rownames(mat))
mat <- mat[idx, ]
rownames(mat) <- mulv_model_genes
mat[is.na(mat)] <- 0

mat_model <- mat_model[rownames(mat), ]


#------ value range ------
cli_alert_info("Value range of query data:")
print(range(mat, na.rm = T))

cli_alert_info("Value range of model data:")
print(range(mat_model, na.rm = T))

#------ Normalize the query matrix ------

if (T) {
    ## In-use
    ## If input query data has >=1 samples.
    ## Multi-mode: rescale for each gene
    mat2 <- sapply(1:nrow(mat), function(i) {
        x <- mat[i, ]
        x <- rescale(x, to = range(mat_model[i, ], na.rm = T))
        x
    })
    mat2 <- t(mat2)
    rownames(mat2) <- rownames(mat)
}
if (T) {
    ## If input query dta only has 1 sample.
    ## find its most similar sample in the model data and then rescale
    ## Single-mode: rescale for each sample
    N_nb <- 1
    dict_most_similar <- sapply(1:ncol(mat), function(j) {
        mat_j <- mat[, j]
        mat_j_nb <- sort(apply(mat_model, 2, function(x) cor(x, mat_j)), decreasing = T) %>% head(N_nb)
        paste0(gtools::mixedsort(names(mat_j_nb)), collapse = ",")
    })
    names(dict_most_similar) <- colnames(mat)
    enframe(dict_most_similar, name = "QuerySample", value = "ModelSample") %>%
        arrange(ModelSample) %>%
        write_csv(file.path(dir_res, "dict.most_similar_sample_in_model.csv"))

    mat3 <- sapply(1:ncol(mat), function(j) {
        mat_j <- mat[, j]
        mat_j_nb <- sort(apply(mat_model, 2, function(x) cor(x, mat_j)), decreasing = T) %>% head(N_nb)
        dict_most_similar <- c(dict_most_similar, head(names(mat_j_nb), 1))
        mat_j_model <- mat_model[, names(mat_j_nb), drop = F]
        mat_j <- rescale(mat_j, to = range(rowMeans(mat_j_model), na.rm = T))
        mat_j
    })
    colnames(mat3) <- colnames(mat)
}

if (T) {
    ## Single-mode: rescale for each sample
    mat_model_consensus <- rowMeans(mat_model, na.rm = T)
    mat4 <- sapply(1:ncol(mat), function(j) {
        mat_j <- mat[, j]
        mat_j <- rescale(mat_j, to = range(mat_model_consensus, na.rm = T))
        mat_j
    })
    colnames(mat4) <- colnames(mat)
}

pal_datatype <- c(
    "Query" = "grey",
    "Query_MultiNorm" = "blue",
    "Query_SingleNorm" = "limegreen",
    "Query_SingleNormCss" = "orange",
    "Model" = "black"
)
df_expr_diag <- rbind(
    as.data.frame(mat) %>% rownames_to_column("Gene") %>% gather("Sample", "Value", -Gene) %>% mutate(data = "Query"),
    as.data.frame(mat2) %>% rownames_to_column("Gene") %>% gather("Sample", "Value", -Gene) %>% mutate(data = "Query_MultiNorm"),
    as.data.frame(mat3) %>% rownames_to_column("Gene") %>% gather("Sample", "Value", -Gene) %>% mutate(data = "Query_SingleNorm"),
    as.data.frame(mat4) %>% rownames_to_column("Gene") %>% gather("Sample", "Value", -Gene) %>% mutate(data = "Query_SingleNormCss"),
    as.data.frame(mat_model) %>% rownames_to_column("Gene") %>% gather("Sample", "Value", -Gene) %>% mutate(data = "Model")
)
df_expr_diag$data <- factor(df_expr_diag$data, levels = c(
    "Query", "Query_MultiNorm",
    "Query_SingleNorm", "Query_SingleNormCss", "Model"
))
p <- ggplot(df_expr_diag, aes(x = Value, y = Gene, color = data)) +
    geom_point(shape = 21, position = position_jitterdodge(jitter.width = .1)) +
    scale_y_discrete(limits = mulv_model_genes) +
    scale_color_manual(values = pal_datatype)
ggsave(
    file.path(dir_res, "predQC.stats.mean_expression.pdf"),
    p + theme(legend.position = "top"),
    width = 5, height = 8, useDingbats = F
)

df_expr_diag_mean <- list(
    enframe(rowMeans(mat, na.rm = T), name = "Gene", value = "Query"),
    enframe(rowMeans(mat_model, na.rm = T), name = "Gene", value = "Model"),
    enframe(rowMeans(mat2, na.rm = T), name = "Gene", value = "Query_MultiNorm"),
    enframe(rowMeans(mat3, na.rm = T), name = "Gene", value = "Query_SingleNorm"),
    enframe(rowMeans(mat4, na.rm = T), name = "Gene", value = "Query_SingleNormCss")
) %>% Reduce(full_join, .)
write_csv(
    df_expr_diag_mean,
    file.path(dir_res, "predQC.stats.mean_expression.csv")
)

source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
library(ggrepel)
for (datatype in c("Query", "Query_MultiNorm", "Query_SingleNorm", "Query_SingleNormCss")) {
    cor_test_combo_obj <- cor_test_combo(df_expr_diag_mean$Model, df_expr_diag_mean[[datatype]])
    p <- df_expr_diag_mean %>%
        ggplot(., aes_string(x = "Model", y = datatype)) +
        geom_point(shape = 21, size = 2)
    p <- p + labs(title = report_cor_test_combo(cor_test_combo_obj)) +
        theme(aspect.ratio = 1) +
        geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
        geom_smooth(method = "lm")
    p <- p + geom_text_repel(aes(label = Gene),
        hjust = -0.2, vjust = 0,
        min.segment.length = unit(0, "lines")
    )
    ggsave(
        file.path(dir_res, sprintf("predQC.scatterplot.stats.mean_expression.%s.pdf", datatype)),
        p,
        width = 5, height = 5, useDingbats = F
    )
    write.csv(
        cor_test_combo_obj,
        file.path(dir_res, sprintf("predQC.correlation.stats.mean_expression.%s.csv", datatype))
    )
}
#------ predict the risk score ------
message("Predicting the risk score...")
for (datatype in c("Query", "Query_MultiNorm", "Query_SingleNormCss", "Query_SingleNorm")) {
    mat_use <- switch(datatype,
        "Query" = mat,
        "Query_MultiNorm" = mat2,
        "Query_SingleNorm" = mat3,
        "Query_SingleNormCss" = mat4
    )
    pred_prob <- predict(mulv_model, as.data.frame(t(mat_use)), type = "response")
    stopifnot(identical(rownames(df), names(pred_prob)))
    df[[sprintf("risk_score_%s", datatype)]] <- pred_prob
    df[[sprintf("risck_score_zscore_%s", datatype)]] <- as.numeric(scale(pred_prob))
}
df$risk_score <- df$risk_score_Query_MultiNorm
write_csv(df, file.path(dir_res, "df_with_risk_score.csv"))


#------ - HR - KM plot ------
suggested_cutoff <- 0.5
for (datatype in c("Query", "Query_MultiNorm", "Query_SingleNorm", "Query_SingleNormCss")) {
    riskscore_str <- sprintf("risk_score_%s", datatype)
    cli_h2(riskscore_str)

    for (i in seq_along(eval_x_choices)) {
        # for (i in 2){
        status_str <- eval_y_choices[i]
        time_str <- eval_x_choices[i]
        cat("\n-- Examine", status_str, "with time", time_str, "\n")

        # j: mimic the dead/alive
        # fixed=0.5 is the most making sense.
        # for (j in c("fixed", "median", "suggested", "mean")) {
        for (j in c("median", "mean")) {
            df_viz <- df
            # df_viz$riskscore <- scale(riskscore)
            df_viz$riskscore <- df_viz[[riskscore_str]]
            df_viz$time <- df_viz[[time_str]]

            prob_cutoff <- switch(j,
                fixed = 0.5,
                median = median(df_viz$riskscore),
                mean = mean(df_viz$riskscore),
                suggested = suggested_cutoff
            )

            df_viz$status <- df_viz[[status_str]]
            print(table(df_viz$status))
            df_viz <- df_viz %>% dplyr::filter(status %in% c(0, 1))
            p <- ggplot(df_viz, aes(x = factor(status), y = riskscore)) +
                geom_boxplot(outlier.shape = NA) +
                geom_quasirandom() +
                stat_compare_means() +
                theme_pubr()
            # print(p)
            # wilcox.test(riskscore ~ df_viz$status)
            tapply(df_viz$riskscore, df_viz$status, mean)

            pdf(
                file.path(dir_res, sprintf(
                    "eval_boxplot.%s.%s.pdf", status_str, datatype
                )),
                width = 6, height = 6, useDingbats = F
            )
            print(p)
            dev.off()

            res <- coxph(Surv(time, status) ~ riskscore, data = df_viz)
            # print(summary(res))
            parse_univariate_coxph <- function(x) {
                x <- summary(x)
                p.value <- signif(x$wald["pvalue"], digits = 2)
                wald.test <- signif(x$wald["test"], digits = 2)
                beta <- signif(x$coef[1], digits = 2) # coeficient beta
                exp_beta <- signif(exp(x$coef[1]), digits = 2)
                HR <- signif(x$coef[2], digits = 2) # exp(beta)
                HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
                HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
                HR <- paste0(
                    HR, " (",
                    HR.confint.lower, "-", HR.confint.upper, ")"
                )
                res <- c(beta, exp_beta, HR, wald.test, p.value)
                names(res) <- c(
                    "beta", "exp_beta", "HR (95% CI for HR)", "wald.test",
                    "p.value"
                )
                return(res)
            }
            coxph_report <- parse_univariate_coxph(res)
            print(coxph_report["p.value"])
            p2 <- try(survminer::ggforest(res, data = df_viz))
            # print(p2)

            # as.numeric(df_viz$riskscore > prob_cutoff)
            # df_viz$grouping <- ifelse(df_viz$riskscore>median(df_viz$riskscore), 'high', 'low')
            df_viz$grouping <- ifelse(df_viz$riskscore > prob_cutoff, "high", "low")
            fit <- survfit(Surv(time, status) ~ grouping, data = df_viz)
            pval <- surv_pvalue(fit)$pval
            print(pval)
            ggsurvplot_p <- ggsurvplot(fit,
                pval = TRUE,
                # conf.int = T,
                risk.table = T, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                # linetype = "strata", # Change line type by groups
                # surv.median.line = "hv", # Specify median survival
                # censor.shape = 124,
                color = "grouping",
                palette = c("high" = "#811C9A", "low" = "#53D43F"),
                ggtheme = theme_pubr()
            )

            ggsurvplot_p$plot <- ggsurvplot_p$plot +
                labs(
                    x = "days", y = sprintf("%s probability", status_str),
                    caption = sprintf("HR=%s coxph P=%s", coxph_report["HR (95% CI for HR)"], coxph_report["p.value"])
                )
            library(patchwork)
            print(ggsurvplot_p, newpage = F)
            pdf(
                file.path(dir_res, sprintf(
                    "eval_survplot.%s.group_by_%s.%s.pdf", status_str, j, datatype
                )),
                width = 6, height = 6, useDingbats = F
            )
            print(ggsurvplot_p, newpage = F)
            dev.off()
            write_csv(
                df_viz[, c("time", "status", "grouping")],
                file.path(
                    dir_res, sprintf("DATAFRAME.%s.group_by_%s.%s.csv", status_str, j, datatype)
                )
            )
            # print(ggsurvplot_p$plot + p)
            rm(ggsurvplot_p)
            rm(fit)
            rm(p)
        }
    }
}


cat("[done]")
timestamp()
