suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Apply the 13-gene model to other cohorts with response data binary outcome
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: text
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
suppressPackageStartupMessages({
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
    library(pROC)
})
#------ parse command line arguments ------
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    study_name <- cmdargs[1]
    f_type <- cmdargs[2]
} else {

    study_name <- "METABRIC"
    f_type <- "CHEMOyes"
    study_name <- "SCANB"
    f_type <- "dUTP"
}

cat("[done]")
timestamp()
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
    sprintf("%s.multivariable_logistic_regression_model.rds", "pCR_status")
))
mulv_model_genes <- read_lines(file.path(
    dir_model,
    sprintf("%s.multivariable_logistic_regression_model_genes.txt", "pCR_status")
))
mulv_model_genes <- intersect(names(mulv_model$coefficients), mulv_model_genes)
#------ model gene matrix ------
mat_model <- read_rds(file.path("/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/psbulk/", "normalized.matrix.rds"))
#------ prepare the input expression data ------
response_str <- "response"
response_lvs <- c("R", "NR")

f_mat <- file.path(
    "/volumes/USR1/yyan/shared/data",
    study_name, f_type,
    "mat.rds"
)
f_df <- file.path(
    "/volumes/USR1/yyan/shared/data",
    study_name, f_type,
    "clinical.df.rds"
)

if (study_name == "2021AyseBassez") {
    study_name <- "2021AyseBassez"
    f_type <- "TNBC_preTx"
    f_mat <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("psbulk_normalized_matrix.%s.rds", f_type)
    )
    f_df <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("df_samplemeta_psbulk.%s.rds", f_type)
    )
    response_str <- "response"
    response_lvs <- c("R", "NR")
}
if (study_name == "artemis_pretx") {
    ## apply to its own
    study_name <- "artemis_pretx"
    f_type <- "psbulk_allcells"
    f_mat <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/psbulk/normalized.matrix.rds"
    f_df <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/psbulk/patient_metainfo.rds"
    response_str <- "pCR_status"
    response_lvs <- c("pCR", "RD")
}
if (study_name == "2024Shiao") {
    study_name <- "2024Shiao"
    f_type <- "TNBC_baseline_cancercells"
    f_mat <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("psbulk_normalized_matrix.%s.rds", f_type)
    )
    f_df <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("df_samplemeta_psbulk.%s.rds", f_type)
    )

}

if (study_name %in% c("BrighTNess")) {
    response_str <- "Response"
    response_lvs <- c("pCR", "RD")
}

if (study_name %in% c("NCPark2020")) {
    response_str <- "Response"
    response_lvs <- c("pCR", "RD")

}

if (study_name == "ISPY990") {
    study_name <- "ISPY990"
    response_str <- "Response"
    response_lvs <- c("pCR", "RD")


    f_mat <- file.path(
        "/volumes/USR1/yyan/shared/data",
        study_name, f_type,
        "mat.rds"
    )
    f_df <- file.path(
        "/volumes/USR1/yyan/shared/data",
        study_name, f_type,
        "clinical.df.rds"
    )
}


# # Only survival data; no response data
if (study_name == "METABRIC") {
    # metabric's rest dataset
    f_mat <- file.path("/volumes/USR1/yyan/shared/data/METABRIC/ruli/", "tnbc_chemoYES_expression.mat.rds")
    f_df <- file.path("/volumes/USR1/yyan/shared/data/METABRIC/ruli/", "tnbc_chemoYES_survival.df.rds")
}


if (study_name == "SCANB") {
    # SCANB TNBC updated 2025-03-04
    dir_lib <- file.path("/volumes/USR1/yyan/shared/data/SCANB_TNBC_chemo/stringtie_fpkm_gene_data_unadjusted/", f_type)
    f_df <- file.path(dir_lib, "scanb_tnbc.df_clinical.rds")
    f_mat <- file.path(dir_lib, "scanb_tnbc.expr_mat.rds")
}

dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas/summary_fig/ML13gene_app",
    study_name, "ml_13gene", f_type
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
    print(range(mat)) # 4.682002 14.638106
    df$response <- ifelse(df$OS_STATUS != 0, "NR", "R")
}

if (study_name == "SCANB") {
    df <- as.data.frame(df)
    rownames(df) <- df$GEX.assay
    eval_x_choices <- c("DRFi_days", "OS_days", "RFi_days", "BCFi_days")
    eval_y_choices <- c("DRFi_event", "OS_event", "RFi_event", "BCFi_event")
    df$response <- ifelse(df$OS_event != 0, "NR", "R")
}


if (study_name == "artemis_pretx") {
    rownames(df) <- df$patient
    df <- df[colnames(mat), ]
}

if (study_name == "ISPY990") {
    df <- as.data.frame(df)
    rownames(df) <- df[["Patient Identifier"]]
    df$response <- df$Response
    mat[is.na(mat)] <- 0
    # mat <- log2(mat + 1)
}

if (study_name == "NCPark2020") {
    df <- as.data.frame(df)
    rownames(df) <- df$sample_id
}

if (study_name == "BrighTNess") {
    df <- as.data.frame(df)
    rownames(df) <- df$patient_id
}
df$response <- df[[response_str]]
df <- df %>% filter(response %in% response_lvs)
df$response <- factor(df$response, levels = response_lvs)
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
    ## In use
    ## If input query data has >=1 samples
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
    width = 5, height = 7, useDingbats = F
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
for (datatype in c("Query", "Query_MultiNorm", "Query_SingleNorm", "Query_SingleNormCss")) {
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
#------ compare ------
write_csv(df, file.path(dir_res, "df_with_risk_score.csv"))
for (datatype in c("Query", "Query_MultiNorm", "Query_SingleNorm", "Query_SingleNormCss")) {
    p <- df %>%
        ggplot(., aes_string(x = "response", y = sprintf("risk_score_%s", datatype))) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom() +
        stat_compare_means()
    ggsave(file.path(dir_res, sprintf("boxplot.risk_score_vs_response_%s.pdf", datatype)), p,
        width = 2.5, height = 3, useDingbats = F
    )

}
for (datatype in c("Query", "Query_MultiNorm", "Query_SingleNorm", "Query_SingleNormCss")) {
    roc_obj <- roc(response = df$response, predictor = df[[sprintf("risk_score_%s", datatype)]])
    pdf(file.path(dir_res, sprintf("ROC.risk_score_vs_response_%s.pdf", datatype)), useDingbats = F)
    plot.roc(
        roc_obj,
        main = sprintf("AUC=%.3f", roc_obj$auc), xlim = c(1, 0), ylim = c(0, 1), asp = 1
    )
    dev.off()
    print(roc_obj$auc)
}


if (study_name == "2021AyseBassez") {
    df$response <- case_when(df$response == "n/a" ~ NA,
        .default = as.character(df$response)
    )
    cohort_opts <- unique(df$cohort)
    df$logic_has_cancer_cells <- ifelse(df$sample_has_cancer_cells > 0, "yes", "no")
    p_list <- lapply(cohort_opts, function(coh) {
        df_s <- df %>%
            dplyr::filter(cohort == coh)
        ggplot(df_s, aes(x = response, y = risk_score)) +
            geom_boxplot(outlier.shape = NA) +
            stat_compare_means(comparisons = list(c("R", "NR"))) +
            geom_quasirandom(aes(color = logic_has_cancer_cells)) +
            scale_x_discrete(limits = c("R", "NR")) +
            labs(title = coh, color = "with cancer cells") +
            scale_color_manual(values = c("yes" = "#0099F9", "no" = "orange4")) +
            coord_cartesian(ylim = c(0, 1))
    })

    p <- patchwork::wrap_plots(p_list, nrow = 1, guides = "collect")
    p
    ggsave(file.path(dir_res, "boxplot.risk_score_vs_response.pdf"), p,
        width = 5, height = 3, useDingbats = F
    )
}

if (study_name == "artemis_pretx") {
    # "pCR_status" "RCB_status" "archetype"

    pdf(file.path(dir_res, "boxplot.risk_score_vs_pCR_status.pdf"),
        width = 2.5, height = 3, useDingbats = F
    )
    (ggplot(df, aes_string(x = "pCR_status", y = "risk_score")) +
        geom_boxplot(outlier.shape = NA) +
        stat_compare_means(comparisons = list(c("pCR", "RD"))) +
        geom_quasirandom()) %>% print()
    dev.off()

    pdf(file.path(dir_res, "boxplot.risk_score_vs_RCB_status.pdf"),
        width = 3.5, height = 3, useDingbats = F
    )
    (ggplot(df, aes_string(x = "RCB_status", y = "risk_score")) +
        geom_boxplot(outlier.shape = NA) +
        stat_compare_means() +
        geom_quasirandom()) %>% print()
    dev.off()

    pdf(file.path(dir_res, "boxplot.risk_score_vs_archetype.pdf"),
        width = 3, height = 3, useDingbats = F
    )
    (ggplot(df, aes_string(x = "archetype", y = "risk_score")) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom() +
        scale_x_discrete(limits = paste0("ARC", 1:4)) +
        stat_compare_means()) %>% print()
    dev.off()
}

cat("[done]")
timestamp()
