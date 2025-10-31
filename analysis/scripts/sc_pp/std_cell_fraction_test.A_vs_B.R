#---------------------------
# Perform cell fraction test at two conditions (e.g., pCR vs RD) ----
#---------------------------
library(Seurat)
library(tidyverse)
library(readr)
library(ggbeeswarm)
suppressPackageStartupMessages({
    library(ggpubr)
    library(forcats)
    library(gtools)
    library(scales)
    library(ggbeeswarm)
    library(cli)
    library(glue)
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
})
theme_set(theme_pubr(legend = "right"))
N_CELLS_MORETHAN <- 100

cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    f_in <- cmdargs[1]
    N_CELLS_MORETHAN <- as.numeric(cmdargs[2])
    on_what <- cmdargs[3]
    by_what <- cmdargs[4]
    patient_str <- cmdargs[5]
    a <- cmdargs[6]
    b <- cmdargs[7]

} else {
}

if (T) {
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    #------ pallete ------
    pal_celltypes <- c(
        "Tumor" = "#f90c71",
        "Epi_Basal" = "slategrey",
        "Epi_LumSec" = "burlywood4",
        "Epi_LumHR" = "antiquewhite4",
        "Mye" = "#009E73",
        "T" = "#AB51E3",
        "B" = "#56B4E9",
        "Fibro" = "#E69F00",
        "Endo" = "#D55E00",
        "Peri" = "#F0E442"
    )
    pal_celltypes2 <- c(
        "Tumor" = "#f90c71",
        "Epithelial" = "#999999",
        "Mye" = "#009E73",
        "T" = "#AB51E3",
        "B" = "#56B4E9",
        "Fibro" = "#E69F00",
        "Endo" = "#D55E00",
        "Peri" = "#F0E442"
    )
    pal_pcr <- c(
        "pCR" = "#53D43F",
        "RD" = "#811C9A",
        "Unknown" = "black",
        "Excluded" = "grey",
        "Removed" = "ghostwhite"
    )
    pal_study <- c(
        "normal" = "#219ebc",
        "tnbc" = "#ffb703"
    )
}

#----------- START -----------
try(print(table(df_cellmeta[[on_what]], useNA = "always")))
print(class(df_cellmeta[[on_what]]))

# try(print(table(fct_drop(df_cellmeta[[on_what]]), useNA = "always")))

stopifnot(all(c(on_what, by_what, patient_str) %in% colnames(df_cellmeta)))

if (on_what == "celltype") {
    df_cellmeta$celltype <- ruok::replace_vector(as.character(df_cellmeta$celltype), c("Epi_Basal" = "Epithelial", "Epi_LumHR" = "Epithelial", "Epi_LumSec" = "Epithelial"))
}
#------------
# Remove patients with too few cells
#------------
dir_res <- file.path(dirname(f_in), sprintf("cell_fraction_test.%s.%s", on_what, by_what))
dir_res <- file.path(
    dirname(f_in), sprintf("cell_fraction_test.%s.%s", on_what, by_what),
    sprintf("keep_patients_Ncell_gt_%s", N_CELLS_MORETHAN)
)
fs::dir_create(dir_res)

dict_pat2by <- df_cellmeta[, c(patient_str, by_what)] %>%
    unique() %>%
    deframe()
df_cellmeta_use <- df_cellmeta %>% dplyr::filter(.data[[by_what]] %in% c(a, b))
try(print(table(fct_drop(df_cellmeta_use[[on_what]]), useNA = "always")))
# df <- df_cellmeta %>%
#   dplyr::filter(.data[[by_what]] %in% c(a, b) ) %>%
#   dplyr::select_at(c(on_what, patient_str, by_what)) %>%
#   dplyr::count(across(all_of(c(on_what, patient_str))), name = 'n') %>%
#   dplyr::add_count(across(all_of(c(patient_str))), wt = n, name = 'N')
# df <- df %>% dplyr::mutate(frac = n / N)

patient_N <- c(table(df_cellmeta_use[[patient_str]]))
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")

quick_summary_stats(patient_N) %>% write.csv(file.path(dir_res, "Ncell_patient_all.stats.csv"))



tab_n <- table(df_cellmeta_use[[patient_str]], df_cellmeta_use[[on_what]])
df <- tab_n %>%
    as.data.frame.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(patient_str)
df <- pivot_longer(df, cols = !patient_str, names_to = on_what, values_to = "n")
df$N <- patient_N[as.character(df[[patient_str]])]
df$frac <- df$n / df$N
df[[by_what]] <- dict_pat2by[df[[patient_str]]]

str(unique(df[[patient_str]]))
# view(df)
tail(df)

## remove patients with few cells
print(sum(patient_N <= 100))
if (sum(patient_N <= N_CELLS_MORETHAN) > 0) {
    low_patient <- names(patient_N[patient_N <= N_CELLS_MORETHAN])

    write_lines(
        low_patient,
        file.path(dir_res, "remove_patients_with_low_cell_counts.txt")
    )
    cli::cli_alert_info(sprintf("Removing patients with low cell counts: %s", length(low_patient)))
    df <- df %>%
        dplyr::filter(!(.data[[patient_str]] %in% low_patient))
}

if (nrow(df) == 0) {
    stop(sprintf("None of the patients have >%s cells", N_CELLS_MORETHAN))
}
#-------------------------- Viz --------------------------
cli_h1('Viz')

#------ plan how patients are ordered ------
on_what_prior <- NULL
pal_use <- init_pal_d(df[[on_what]], pal = "circus")

if (on_what == "celltype") {
    df$celltype <- factor(df$celltype, levels = names(pal_celltypes2))
    df$celltype <- forcats::fct_drop(df$celltype)
    on_what_prior <- "Tumor"
    pal_use <- pal_celltypes2
}
if (on_what == "cell_state_paper") {
    std_id_lvs <- read_rds("/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/idents_levels_cellstates.rds")
    if (length(setdiff(as.character(unique(df_cellmeta_use$cell_state_paper)), std_id_lvs)) > 0) {
        on_what_lvs <- c(
            intersect(std_id_lvs, as.character(unique(df_cellmeta_use$cell_state_paper))),
            setdiff(as.character(unique(df_cellmeta_use$cell_state_paper)), std_id_lvs)
        )
    } else {
        on_what_lvs <- intersect(std_id_lvs, as.character(unique(df_cellmeta_use$cell_state_paper)))
    }
    
    df[[on_what]] <- factor(
        as.character(df[[on_what]]),
        levels = on_what_lvs
    )
}
if (on_what == "MetaNiche") {
    df[[on_what]] <- factor(
        as.character(df[[on_what]]),
        levels = names(pal_metaniche)
    )
    pal_use <- pal_metaniche
}

if (! "factor" %in% class(df[[on_what]])) {
    message(sprintf("Setting %s as factor with mixedsort levels", on_what))
    df[[on_what]] <- factor(
        as.character(df[[on_what]]),
        levels = gtools::mixedsort(unique(as.character(df[[on_what]])))
    )
}
print(table(df[[on_what]], useNA = "always"))

n_on_what <- length(levels(df[[on_what]]))

if (is.null(on_what_prior)) {
    # patient_order <- gtools::mixedsort( unique(df[[patient_str]]) )
    # names(patient_order) <- patient_order
    on_what_prior <- df %>%
        dplyr::group_by_at(on_what) %>%
        dplyr::summarise(avg_frac = mean(frac)) %>%
        dplyr::top_n(n = 1, wt = avg_frac) %>%
        dplyr::pull(on_what) %>%
        as.character()
}
patient_order <- df %>%
    dplyr::filter(.data[[on_what]] == on_what_prior) %>%
    dplyr::select_at(c(patient_str, "frac")) %>%
    deframe() %>%
    sort()
df[[patient_str]] <- factor(as.character(df[[patient_str]]),
    levels = names(patient_order)
)

#------ barplot of patients A vs B ------
cli_h2("Barplot of patients A vs B")
# stack/fill
for (ggplot_bar_type in c("fill", "stack")) {
    p1 <- ggplot(df, aes_string(x = patient_str, y = "n", fill = on_what)) +
        geom_bar(position = ggplot_bar_type, stat = "identity") +
        facet_wrap(as.formula(sprintf(".~%s", by_what)), scales = "free_x") +
        labs(y = switch(ggplot_bar_type,
            fill = "frac",
            stack = "num"
        )) +
        scale_fill_manual(values = pal_use) +
        rotate_x_text(90) +
        theme(axis.text.x = element_text(size = 4))
    ggsave(
        file.path(dir_res, sprintf(
            "barplot.%s.pdf",
            switch(ggplot_bar_type,
                fill = "frac",
                stack = "num"
            )
        )),
        p1,
        width = 10, height = 5, useDingbats = F
    )
}

#------ barplot summarized A v.s. B ------
cli_h2("Barplot summarized A vs B")
head(df)
df_summary <- df %>%
    dplyr::group_by_at(c(by_what, on_what)) %>%
    dplyr::summarise(value = mean(n))
head(df_summary)
write_csv(df_summary, file.path(dir_res, "average.cell_num.csv"))
lapply(split(df_summary, df_summary[[by_what]]), function(x) {
    x$frac <- x$value / sum(x$value)
    x
}) %>% do.call("rbind", .) -> df_summary
write_csv(df_summary, file.path(dir_res, "summary.pseudo_profile.csv"))

# fill/stack
for (ggplot_bar_type in c("fill", "stack")) {
    p <- ggplot(df_summary, aes_string(x = by_what, y = "value", fill = on_what)) +
        geom_bar(position = ggplot_bar_type, stat = "identity") +
        labs(y = switch(ggplot_bar_type,
            fill = "avg frac",
            stack = "avg num"
        )) +
        scale_fill_manual(values = pal_use)
    ggsave(
        file.path(dir_res, sprintf(
            "barplot.summaried.%s.pdf",
            switch(ggplot_bar_type,
                fill = "frac",
                stack = "num"
            )
        )),
        p,
        width = 4, height = 5, useDingbats = F
    )
}
## lolipop plot
df_summary$on_what <- df_summary[[on_what]]
df_summary$by_what <- df_summary[[by_what]]
p_df <- df_summary %>% 
    dplyr::group_by_at('on_what') %>%
    dplyr::mutate(diff_fraction = frac[by_what == b] - frac[by_what == a]) %>%
    dplyr::select(on_what, diff_fraction) %>%
    unique()
p_df$color <- ifelse(p_df$diff_fraction >= 0, b, a)

p <- ggplot(p_df, aes(
    y = reorder(on_what, diff_fraction),
    x = diff_fraction
)) +
    geom_col(aes(fill = color)) +
    scale_x_continuous(
        labels = scales::percent_format(accuracy = 1),
        breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_fill_manual(values = pal_pcr) +
    rremove("legend") +
    theme(legend.position = "right")
# p
ggsave(
    file.path(
        dir_res,
        sprintf("summary.pseudo_profiles.lolipop.cell_fraction_diff.PCR.colored.pdf")
    ),
    p,
    width = 5,
    height = length(levels(p_df[['on_what']])) * 0.4,
    useDingbats = F, limitsize = F
)
df_summary$on_what  <- NULL
df_summary$by_what  <- NULL


#------------------ ~~~ barplot real frac A v.s. B ~~~ --------------------
cli_h1("barplot real frac A v.s. B")

df_summary_frac <- df %>%
    dplyr::group_by_at(c(by_what, on_what)) %>%
    dplyr::summarise(value = mean(frac))
write_csv(df_summary_frac, file.path(dir_res, "average.cell_frac.csv"))
p <- ggplot(df_summary_frac, aes_string(x = by_what, y = "value", fill = on_what)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(y = "avg frac") +
    scale_fill_manual(values = pal_use)
ggsave(
    file.path(dir_res, sprintf(
        "barplot.summaried.%s.pdf",
        "real_frac"
    )),
    p,
    width = 4, height = 5, useDingbats = F
)

#------ heatmap summarized A vs B ------
cli_h2("Heatmap summarized A vs B")
mat_summary <- df_summary %>%
    dplyr::select(all_of(c(on_what, by_what, "value"))) %>%
    pivot_wider(names_from = on_what, values_from = "value") %>%
    column_to_rownames(by_what) %>%
    as.matrix()

chi_res <- chisq.test(as.table(mat_summary))
chi_res$residuals %>%
    as.data.frame() %>%
    ggplot(aes(x = Freq, y = Var1)) +
    geom_col(aes(fill = Freq >= 0)) +
    facet_wrap(~Var2, nrow = 2) +
    scale_fill_manual(values = c(`TRUE` = "orange", `FALSE` = "grey")) +
    ggthemes::theme_base() +
    rremove("legend") +
    geom_vline(xintercept = 0, lty = "dashed") +
    labs(
        x = "pearson residuals", y = by_what,
        caption = sprintf(
            "chi-square (%s, %s)=%s pval=%s",
            chi_res$parameter, sum(chi_res$observed), round(chi_res$statistic, 2), scientific(chi_res$p.value, digits = 2)
        )
    )

#------ wilcox test ------
cli_h2("Wilcox test")
head(df)

head(df)
tmp <- split(df, df[[on_what]])
manual_test_res <- lapply(tmp, function(xx) {
    # cat(dim(xx), '..\n')
    ggpubr::compare_means(data = xx, formula = as.formula(sprintf("%s~%s", "frac", by_what)))
})
for (x in names(manual_test_res)) {
    manual_test_res[[x]][[on_what]] <- x
}
rm(x)
manual_test_res <- do.call("rbind", manual_test_res)
manual_test_res$p.adj <- p.adjust(manual_test_res$p, method = "fdr")
manual_test_res$q.format <- paste0("q=", as.character(signif(manual_test_res$p.adj, digits = 3)))
manual_test_res$p.format <- paste0("p=", as.character(signif(manual_test_res$p, digits = 3)))

manual_test_res[[on_what]] <- factor(as.character(manual_test_res[[on_what]]), levels = levels(df[[on_what]]))
library(ggpubr)
library(rstatix)

manual_qval_text <- ruok::pretty_table2str(
    deframe(manual_test_res[, c(on_what, "p.adj")]) %>% signif(., digits = 3)
)
manual_qval <- deframe(manual_test_res[, c(on_what, "p.adj")])


# pal_use <- pal_study
pal_use <- pal_pcr
p <- ggplot(df, aes_string(x = on_what, y = "frac", color = by_what)) +
    geom_quasirandom(
        pch = 21, color = "black", aes_string(fill = by_what), size = 2,
        stroke = 0.5,
        varwidth = TRUE, dodge.width = 0.9
    ) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.9), fill=NA) +
    # geom_jitter(pch=21, color='black', aes_string(fill=by_what), size=1) +
    scale_color_manual(values = pal_use) +
    scale_fill_manual(values = pal_use) +
    facet_wrap(as.formula(sprintf(".~%s", on_what)),
        nrow = 1,
        scales = "free",
        labeller = as_labeller(manual_qval_text)
    ) +
    rotate_x_text(45)
# theme(strip.text = element_text(size = 4), axis.text = element_text(size=4))
# p
ggsave(file.path(dir_res, "boxplot.test.faceted.pdf"), p,
    width = n_on_what * 1.4 + 0.5, # ceiling(sqrt(n_on_what)) * 1.5,
    height = 5, # ceiling(sqrt(n_on_what)) * 1.5,
    useDingbats = F, limitsize = F
)

p <- ggplot(df, aes_string(x = on_what, y = "frac")) +
    geom_quasirandom(
        pch = 21, color = "black", aes_string(fill = by_what), size = 2,
        stroke = 0.5,
        varwidth = TRUE, dodge.width = 0.9
    ) +
    geom_boxplot(
        outlier.shape = NA, aes_string(color = by_what),
        position = position_dodge(0.9), fill=NA
    ) +
    scale_color_manual(values = pal_use) +
    scale_fill_manual(values = pal_use) +
    stat_pvalue_manual(
        manual_test_res,
        label = after_stat("p.format"),
        x = on_what, bracket.nudge.y = 1, y.position = 1, remove.bracket = F
    ) +
    stat_pvalue_manual(
        manual_test_res,
        label = after_stat("q.format"),
        x = on_what, bracket.nudge.y = 1, y.position = 1.1, remove.bracket = F
    ) +
    rotate_x_text(45)
# p
ggsave(file.path(dir_res, "boxplot.test.compact.pdf"), p,
    width = n_on_what * 1.1 + 0.5,
    height = 5, useDingbats = F, limitsize = F
)
p <- ggplot(df, aes_string(x = on_what, y = "frac", fill = by_what)) +
    geom_quasirandom(
        pch = 21, color = "black", 
        # aes_string(fill = by_what), 
        size = 2,
        stroke = 0.5,
        varwidth = TRUE, dodge.width = 0.9
    ) +
    ## add mean points
    stat_summary(
        fun = mean,
        geom = "point", shape = 18, size = 3,
        # aes_string(fill = by_what), 
        color = "orange",
        position = position_dodge(0.9)
    ) +
    stat_summary(
        fun = mean, geom = "errorbar",
        aes(ymax = after_stat(y), ymin = after_stat(y)), 
        position = position_dodge(width = 0.9)
    ) + scale_fill_manual(values = pal_use) + scale_color_manual(values = pal_use) + 
    rremove("legend") 
    # + 
    # stat_pvalue_manual(
    #     manual_test_res,
    #     label = after_stat("p.format"),
    #     x = on_what, bracket.nudge.y = 1, y.position = 1, remove.bracket = F
    # )
    #  +
    # stat_pvalue_manual(
    #     manual_test_res,
    #     label = after_stat("q.format"),
    #     x = on_what, bracket.nudge.y = 1, y.position = 1.1, remove.bracket = F
    # )
p
ggsave(file.path(dir_res, "quasirandom.test.compact.pdf"), p,
    width = n_on_what * 1.1 + 0.5,
    height = 5, useDingbats = F, limitsize = F
)


#------ Check other test results ------
cli_h2("Check other test results")
for (difftest_method in c("wilcoxon", "dunn", "ttest", "kruskal")) {
    test_use <- switch(difftest_method,
        wilcoxon = wilcox_test,
        dunn = dunn_test,
        ttest = t_test,
        kruskal = kruskal_test,
        signtest = sign_test,
        anovatest = anova_test
    )
    for (pvaladj in c("fdr", "holm", "hochberg", "hommel")) {
        glue("Test: {difftest_method} with pval adjustment: {pvaladj}")
        stat_test_df <- df %>%
            group_by_at(on_what) %>%
            test_use(
                formula = as.formula(sprintf("%s ~ %s", "frac", by_what))
            ) %>%
            adjust_pvalue(method = pvaladj) %>%
            add_significance("p.adj") %>%
            as.data.frame() %>%
            dplyr::select(on_what, p, p.adj, p.adj.signif) %>%
            mutate(padjmethod = pvaladj, diffmethod = difftest_method)

        write_csv(stat_test_df, file.path(dir_res, sprintf(
            "report.cell_fraction_test.%s.%s.%s.%s.csv",
            on_what, by_what, difftest_method, pvaladj
        )))
        write_rds(stat_test_df, file.path(dir_res, sprintf(
            "report.cell_fraction_test.%s.%s.%s.%s.rds",
            on_what, by_what, difftest_method, pvaladj
        )))
    }
}






#------ Lolipop plot ------
cli_h2("Lolipop plot")

df$on_what <- df[[on_what]]
df$by_what <- df[[by_what]]
mm_viz <- unique(df[[on_what]])
p_df <- df %>%
    dplyr::filter(by_what %in% c("pCR", "RD")) %>%
    dplyr::group_by_at(c("by_what", "on_what")) %>%
    dplyr::summarise(frac = mean(frac)) %>%
    dplyr::group_by_at("on_what") %>%
    dplyr::mutate(diff_fraction = frac[by_what == "RD"] - frac[by_what == "pCR"]) %>%
    dplyr::select(on_what, diff_fraction) %>%
    unique() %>%
    # dplyr::mutate(color=ifelse(diff_fraction>=0, 'RD', 'pCR')) %>%
    dplyr::mutate(abs_diff_fraction = abs(diff_fraction))
p_df$is_sig <- manual_qval[p_df$on_what]
p_df$is_sig[is.na(p_df$is_sig)] <- 1
p_df$is_sig_cat <- ifelse(p_df$is_sig < 0.05, "sig", "ns")
p_df$is_sig_cat <- factor(p_df$is_sig_cat, levels=c('sig', 'ns'))
p_df$color <- ifelse(p_df$diff_fraction >= 0, "RD", "pCR")
p <- ggplot(p_df, aes(
    y = reorder(on_what, diff_fraction),
    x = diff_fraction
)) +
    geom_col(aes(fill = color)) +
    geom_text(aes(x=0, label = signif(p_df$is_sig, 3)),
              color = "black") +
    # scale_y_discrete(limits=rev(c(mm_viz)))  +
    scale_x_continuous(
        labels = scales::percent_format(accuracy = 1),
        breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_fill_manual(values = pal_pcr) +
    # scale_fill_manual(values = c('sig'='goldenrod1', 'ns'='grey')) +
    facet_grid(rows = "is_sig_cat", scales = "free", space = "free") +
    rremove("legend") +
    theme(legend.position = "right")
# p
ggsave(
    file.path(
        dir_res,
        sprintf("lolipop.cell_fraction_diff.PCR.colored.pdf")
    ),
    p,
    width = 5,
    height = length(mm_viz) * 0.4,
    useDingbats = F, limitsize = F
)



#------ Variation Test ------
cli_h2("Variation Test")
keep_robust_group <- function(x, use_prop = 0.9) {
    if (length(x) <= 2) {
        return(x)
    }
    ## remove outliers, keep the top 90% of the values
    q_lo <- (1 - use_prop) / 2
    q_hi <- 1 - (1 - use_prop) / 2
    x_quant <- quantile(x, probs = c(q_lo, q_hi), na.rm = TRUE)
    return(x[x >= x_quant[1] & x <= x_quant[2]])
}

calc_variation_jeong <- function(ref, alt, do_outlier_robust = F, use_prop = 0.9) {
    if (do_outlier_robust) {
        ref <- keep_robust_group(ref, use_prop = use_prop)
        alt <- keep_robust_group(alt, use_prop = use_prop)
    }
    n_ref <- length(ref)
    n_alt <- length(alt)
    variation_matrix <- as.matrix(cbind(
        i = rep(1:n_alt, each = n_ref),
        j = rep(1:n_ref, times = n_alt)
    ))
    # if (identical(ref, alt)) {
    #     idx <- variation_matrix[, 'i'] != variation_matrix[, 'j']
    #     variation_matrix <- variation_matrix[idx, ]
    # }
    variation_matrix_value <- alt[variation_matrix[, "i"]] / ref[variation_matrix[, "j"]]
    # 0/0=NaN, x/0=Inf
    na_idx <- is.na(variation_matrix_value)
    if (all(na_idx)) {
        return(rep(0, nrow(variation_matrix)))
    }
    variation_matrix_value <- variation_matrix_value[!na_idx]
    # [0, Inf] -> log2 -> [-Inf, Inf]
    variation_matrix_value <- log2(variation_matrix_value) # 0->-Inf; Inf->Inf
    bad_num_idx <- is.na(variation_matrix_value) | !is.finite(variation_matrix_value)
    good_num <- variation_matrix_value[!bad_num_idx]
    if (length(good_num) == 0) {
        return(rep(0, length(variation_matrix_value)))
    }
    variation_matrix_value[variation_matrix_value == -Inf] <- min(good_num) * 0.9
    variation_matrix_value[variation_matrix_value == Inf] <- max(good_num) * 1.1

    return(variation_matrix_value)
}
calc_variation_yun <- function(ref, alt, do_outlier_robust = F, use_prop = 0.9) {
    if (do_outlier_robust) {
        ref <- keep_robust_group(ref, use_prop = use_prop)
        alt <- keep_robust_group(alt, use_prop = use_prop)
    }
    n_ref <- length(ref)
    n_alt <- length(alt)
    variation_matrix <- as.matrix(cbind(
        i = rep(1:n_alt, each = n_ref),
        j = rep(1:n_ref, times = n_alt)
    ))
    # if (identical(ref, alt)) {
    #     idx <- variation_matrix[, 'i'] != variation_matrix[, 'j']
    #     variation_matrix <- variation_matrix[idx, ]
    # }
    variation_matrix_value <- alt[variation_matrix[, "i"]] - ref[variation_matrix[, "j"]]
    return(variation_matrix_value)
}

mm_viz <- unique(df[[on_what]])
df$on_what <- df[[on_what]]
df$by_what <- df[[by_what]]
for (variation_method in c("yy", "jeong", "yy_robust", "jeong_robust")) {
    cli_h3(sprintf("Variation method: %s", variation_method))
    res <- lapply(mm_viz, function(mm) {
        message(mm)
        ref_v <- df %>%
            dplyr::filter(on_what == mm) %>%
            dplyr::filter(by_what == "pCR") %>%
            dplyr::select(frac) %>%
            deframe()
        alt_v <- df %>%
            dplyr::filter(on_what == mm) %>%
            dplyr::filter(by_what == "RD") %>%
            dplyr::select(frac) %>%
            deframe()
        if (variation_method == "jeong") {
            ref_variation <- calc_variation_jeong(ref = ref_v, alt = ref_v)
            alt_variation <- calc_variation_jeong(ref = ref_v, alt = alt_v)
        }
        if (variation_method == "yy") {
            ref_variation <- calc_variation_yun(ref = ref_v, alt = ref_v)
            alt_variation <- calc_variation_yun(ref = ref_v, alt = alt_v)
        }
        if (variation_method == "jeong_robust") {
            ref_variation <- calc_variation_jeong(
                ref = ref_v, alt = ref_v,
                do_outlier_robust = TRUE, use_prop = 0.9
            )
            alt_variation <- calc_variation_jeong(
                ref = ref_v, alt = alt_v,
                do_outlier_robust = TRUE, use_prop = 0.9
            )
        }
        if (variation_method == "yy_robust") {
            ref_variation <- calc_variation_yun(
                ref = ref_v, alt = ref_v,
                do_outlier_robust = TRUE, use_prop = 0.9
            )
            alt_variation <- calc_variation_yun(
                ref = ref_v, alt = alt_v,
                do_outlier_robust = TRUE, use_prop = 0.9
            )
        }

        print(range(ref_variation))
        print(range(alt_variation))

        # if (variation_method == "yy") {
        #     pval <- wilcox.test(x = alt_variation, y = ref_variation)
        # } else {
        #     pval <- wilcox.test(x = alt_variation, y = ref_variation)
        # }
        # pval <- ks.test(x = alt_variation, y = ref_variation, alternative = "two.sided")
        pval <- wilcox.test(alt_variation, ref_variation)
        pval <- pval$p.value
        return(list(pval = pval, alt_variation = alt_variation, ref_variation = ref_variation))
    })
    names(res) <- mm_viz

    ## boxplot
    pdf(file.path(dir_res, sprintf("boxplot.ks_test.variation_combo.%s.by_%s.pdf", variation_method, by_what)),
        width = 2, height = 4, onefile = T, useDingbats = F
    )
    for (mm in mm_viz) {
        tmp <- data.frame(
            value = c(res[[mm]]$ref_variation, res[[mm]]$alt_variation),
            cat = c(
                rep("ref", length(res[[mm]]$ref_variation)),
                rep("alt", length(res[[mm]]$alt_variation))
            )
        )
        tmp$cat <- factor(tmp$cat, levels = c("ref", "alt"))
        tmp_p <- ggplot(tmp, aes(x = cat, y = value, fill = cat)) +
            geom_boxplot(outlier.shape = NA) +
            scale_fill_manual(values = c(ref = "grey", alt = "blue")) +
            labs(title = mm) +
            geom_hline(yintercept = 0, lty = "dashed", color = "black")
        tmp_p <- tmp_p + rremove("legend") + stat_compare_means()
        print(tmp_p)
    }
    rm(tmp)
    rm(tmp_p)
    dev.off()

    p_df <- lapply(mm_viz, function(mm) {
        o <- res[[mm]]
        v <- o$alt_variation
        data.frame(
            on_what = rep(mm, length(v)),
            alt_variation = v
        )
    })
    p_df <- do.call("rbind", p_df)

    p_pval <- sapply(res, function(o) o$pval, simplify = T)
    p_qval <- p.adjust(p_pval, method = "bonferroni") # FDR adjustment
    names(p_qval) <- names(p_pval) <- names(res)
    p_pval_str <- sapply(p_pval, scales::scientific) %>% ruok::pretty_table2str()
    p_qval_str <- sapply(p_qval, scales::scientific) %>% ruok::pretty_table2str()

    p_sig_str <- p_qval_str
    p_sigval <- p_qval
    p_sigval[is.na(p_sigval)] <- 1
    p_df_avg <- tapply(p_df$alt_variation, p_df$on_what, mean)
    p_df_color <- ifelse(p_df_avg > 0, "RD", "pCR")
    p_sigval <- p_sigval[names(p_df_avg)]
    stopifnot(names(p_sigval) == names(p_df_avg)) # ensure order is the same
    p_df_color[p_sigval >= 0.05] <- "n.s."
    p_df$color <- p_df_color[p_df$on_what]
    p_df$on_what <- factor(as.character(p_df$on_what), levels = rev(mm_viz)) # for viz purpose
    p1 <- ggplot(p_df, aes(
        x = alt_variation,
        y = reorder(on_what, -1 * alt_variation, mean),
        fill = color
    )) +
        geom_boxplot(outlier.shape = NA) +
        # geom_violin()+
        # geom_quasirandom(pch=21, color='black', size=2,
        #                  stroke = 0.5,
        #                  varwidth = TRUE, dodge.width=0.9)+
        # stat_summary(fun=median, geom="point", shape=20, color='gold', size=3) +
        stat_summary(fun = mean, geom = "point", shape = 20, color = "gold", size = 3) +
        scale_x_continuous(position = "top", breaks = scales::pretty_breaks(n = 5)) +
        geom_vline(xintercept = 0, lty = "dashed") +
        # scale_y_discrete(limits=rev(mm_viz), labels=p_pval_str) +
        scale_y_discrete(labels = p_sig_str) +
        facet_grid(rows = "color", scales = "free_y", space = "free_y") +
        labs(x = ifelse(variation_method %in% c("yy", "yy_robust"),
            "proportion of cells RD - pCR",
            "log2(proportion of cells RD/pCR)"
        ), subtitle = "all") +
        scale_fill_manual(values = pal_pcr) +
        rremove("legend")
    if (variation_method %in% c("yy", "yy_robust")) {
        p1 <- p1 + scale_x_continuous(
            labels = scales::percent_format(accuracy = 1),
            # limits = c(-0.5, 0.5),
            breaks = scales::pretty_breaks(n = 5),
            position = "top"
        )
    }
    # print(p1)
    ggsave(
        file.path(
            dir_res,
            sprintf("violin.ks_test.variation_%s.by_%s.pdf", variation_method, by_what)
        ),
        p1,
        width = 6, height = 7, useDingbats = F
    )
    p_df_stats <- enframe(p_df_avg)
    p_df_stats$pval <- p_pval[p_df_stats$name]
    p_df_stats$qval <- p_qval[p_df_stats$name]
    write.csv(p_df_stats,
        file = file.path(dir_res, sprintf("report.ks_test.variation_%s.by_%s.csv", variation_method, by_what))
    )
}
df$on_what <- NULL
df$by_what <- NULL

#------ export ------
write_rds(manual_test_res, file.path(dir_res, "manual_test_res.wilcoxon.rds"))
write_csv(manual_test_res, file.path(dir_res, "manual_test_res.wilcoxon.csv"))

write_rds(df, file.path(dir_res, "df.rds"))
write_csv(df, file.path(dir_res, "df.csv"))

write_csv(x = df_summary, file = file.path(dir_res, "df_summary.csv"))
write_rds(x = df_summary, file = file.path(dir_res, "df_summary.rds"))

# rm(list = ls())

df[, c(patient_str, by_what)] %>%
    unique() %>%
    deframe() %>%
    table() %>%
    write.csv(file.path(dir_res, sprintf("N_%s_by_%s.csv", patient_str, by_what)))

cat("Done\n\n")
