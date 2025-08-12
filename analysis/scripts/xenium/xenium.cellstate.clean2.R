suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Proposing and performing QC (clean2, i.e. at cell state level).
#
# Input:
# - Seurat object (clean1)
# - Consensus MapQuery data frame
# Output:
# - Diagnosis goes to the same directory of css.
# - New object goes to clean2 version of the input Seurat object .
# CellLabel: I use the column `css_id` to always save the best result (i.e., clean2's result in the clean2 object).
#
# For now, I run for each sample. But I am compatible with any.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: Nov 13, 2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
suppressPackageStartupMessages({
    library(readr)
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    theme_set(theme_pubr(base_size = 6, legend = "right") %+replace% theme(axis.ticks.length = unit(0.1, "inch")))
    library(cli)
    library(tictoc)
    library(glue)
    library(scales)
    library(tools)
    library(fs)
    library(nanoparquet)
    library(Seurat)
    library(patchwork)
    library(ruok)
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    f_obj <- cmdargs[1]
    f_df <- cmdargs[2]
} else {
    f_obj <- "/volumes/USR1/yyan/project/tnbc_xenium/data/ART277/cleaned1/objects_split_into_celltype/T/pass.xenium.seurat.rds"
    f_df <- "/volumes/USR1/yyan/project/tnbc_xenium/data/ART277/cleaned1/objects_split_into_celltype/T/mapquery_css_cellstate/consensus.dataframe.rds"
}


#------------------ ~~~ Setup Outdir ~~~ --------------------
cli_h1("Setup Outdir")

dir_obj <- str_replace(dirname(f_obj), pattern = "cleaned1", replacement = "cleaned2")
fs::dir_create(dir_obj)

dir_diag <- dirname(f_df)

#------------------ ~~~ Read-in ~~~ --------------------
cli_h1("Read-in")
obj <- read_rds(f_obj)
print(obj)
df_css <- read_rds(f_df)
print(tail(df_css))

stopifnot("css_id" %in% colnames(df_css))

#------------------ ~~~ colors ~~~ --------------------
std_ident_levels <- read_rds(file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate",
    "pat102/atlas",
    "idents_levels_cellstates.rds"
))
ident_levels <- get_levels(df_css$css_id)
print(ident_levels)
pal_ident <- init_pal_d(df_css$css_id)
#------------------ ~~~ Diagnosis ~~~ --------------------
cli_h1("Diagnosis")


#------ nUMI/nFeature/Prob by nAgreements ------
# basic_metric <- "css_prob"
for (basic_metric in c("nCount_Xenium", "nFeature_Xenium", "css_prob")) {
    cat(basic_metric, "... ")
    p <- df_css %>%
        dplyr::mutate(nagreements = as.factor(nagreements)) %>%
        ggplot(., aes_string(x = "nagreements", y = basic_metric)) +
        geom_violin(scale = "area", fill = "lightgrey", color = NA) +
        stat_mean(color = "black", pch = 16) +
        # scale_x_discrete(limits = factor(1:10)) +
        labs(y = sprintf(basic_metric), x = "num agreements")
    p <- p + stat_compare_means(ref.group = "10", label = "p.signif")
    if (basic_metric %in% c("nCount_Xenium", "nFeature_Xenium")) {
        p <- p + scale_y_log10() + annotation_logticks(sides = "l")
    }
    ggsave(
        file.path(
            dir_diag,
            sprintf("violin.%s.all_cells.pdf", basic_metric)
        ),
        p,
        width = 4, height = 2.5, useDingbats = F
    )
}
cat("\n")


#------ nUMI/nFeature/Prob by nAgreements facted by css_id ------
for (basic_metric in c("nCount_Xenium", "nFeature_Xenium", "css_prob")) {
    cat(basic_metric, ">>> ")
    fs::dir_create(
        file.path(dir_diag, sprintf("violin.%s.faceted_by_ident", basic_metric))
    )

    for (id in ident_levels) {
        cat(id, "...")
        p <- df_css %>%
            dplyr::filter(css_id == id) %>%
            dplyr::mutate(nagreements = as.factor(nagreements)) %>%
            ggplot(., aes_string(x = "nagreements", y = basic_metric)) +
            geom_violin(scale = "area", fill = pal_ident[id], color = NA) +
            stat_mean(color = "black", pch = 16) +
            labs(y = sprintf(basic_metric), x = "num agreements")
        p <- p + stat_compare_means(ref.group = "10", label = "p.signif")
        if (basic_metric %in% c("nCount_Xenium", "nFeature_Xenium")) {
            p <- p + scale_y_log10() + annotation_logticks(sides = "l")
        }
        ggsave(
            file.path(
                file.path(dir_diag, sprintf("violin.%s.faceted_by_ident", basic_metric)),
                sprintf("violin.%s.pdf", id)
            ),
            p,
            width = 4, height = 2.5, useDingbats = F
        )
    }
    cat("\n")
}
cat("\n")


#------------------ ~~~ Proposing QC ~~~ --------------------
cli_h1("Proposing QC")
#------ barplot nagreements ------
p <- ruok::qbarplot_table_cat(df_css$nagreements, name_x = "nAgreements") +
    scale_fill_viridis_d()
ggsave(file.path(dir_diag, "barplot.nagreement.pdf"), p, width = 2, height = 3, useDingbats = F)
#------ Function: find cutoff by testing ------
adhoc_find_nagreement_cutoff <- function(nagreements, y, col_sig = "p") {
    ##
    ## The biggest nagreements that is differnt from nagreements=10
    ## - statistical p<0.05
    ## - signal strength: median is either higher or lower than 0.5 * nagreements=10
    ##
    res <- 10

    df_test <- data.frame(x = nagreements, y = y)
    df_test$x <- factor(df_test$x)
    o_test <- try(compare_means(y ~ x, data = df_test, ref.group = "10", p.adjust.method = "BH"))

    # print(o_test)
    ## corner case
    if ("try-error" %in% class(o_test)) {
        cat("error is catched.")
        return(NA)
    }

    strength_comp <- tapply(df_test$y, df_test$x, median) 
    strength_comp_ref <- strength_comp['10'] %>% as.numeric()
    strength_comp <- strength_comp %>% 
        enframe() %>% 
        dplyr::filter(name != '10') %>%
        dplyr::arrange(desc(name))

    o_test$group2 <- as.numeric(o_test$group2)
    o_test <- dplyr::arrange(o_test, desc(group2))


    i <- 1
    while (i <= nrow(o_test)) {
        is_sig_diff <- FALSE
        is_strong_diff <- FALSE

        if (o_test[[col_sig]][i] < 0.05) {
            is_sig_diff <- TRUE
        }
        if (strength_comp$value[i] < .6 * strength_comp_ref | strength_comp$value[i] > 1.1 * strength_comp_ref) {
            is_strong_diff <- TRUE
        }
        if( is_sig_diff & is_strong_diff ){
            break()
        }
        res <- o_test$group2[i]
        i <- i + 1
    }

    return(res)
}

#------ cutoff of all cells ------
basic_metric_opt <- c("nCount_Xenium", "nFeature_Xenium", "css_prob")
global_nagreements_goodmin_opts <- sapply(basic_metric_opt, function(basic_metric) {
    adhoc_find_nagreement_cutoff(df_css[["nagreements"]], df_css[[basic_metric]])
})
print(global_nagreements_goodmin_opts)

#------ cutoff of each cell id ------

id_nagreements_goodmin_opts <- sapply(ident_levels, function(id) {
    cat(id, "... ")
    df_css_id <- df_css %>% dplyr::filter(css_id == id)
    sapply(basic_metric_opt, function(basic_metric) {
        adhoc_find_nagreement_cutoff(df_css_id[["nagreements"]], df_css_id[[basic_metric]])
    })
})
print(id_nagreements_goodmin_opts)
write_rds(global_nagreements_goodmin_opts, file.path(dir_diag, "global_nagreements_goodmin_opts.rds"))
write_rds(id_nagreements_goodmin_opts, file.path(dir_diag, "id_nagreements_goodmin_opts.rds"))
write.csv(global_nagreements_goodmin_opts, file.path(dir_diag, "global_nagreements_goodmin_opts.csv"))
write.csv(id_nagreements_goodmin_opts, file.path(dir_diag, "id_nagreements_goodmin_opts.csv"))

#------ Refine cutoff for each id ------
# nCount_Xenium nFeature_Xenium        css_prob
#               4               4              10
#                 CD4-TN CD4-TFH CD4-TCM CD4-TREG CD4-TIFN CD8-TEM CD8-TEFF
# nCount_Xenium        9       6      10        9       10       5        2
# nFeature_Xenium      9       6      10        9       10       5        2
# css_prob            10       9      10       10       10       9        9
#                 CD8-TRM CD8-TEXH CD8-TIFN GD-T NK-CD16high NK-CD16low T-prolif
# nCount_Xenium         4       NA        3    2           2          8        8
# nFeature_Xenium       5       NA        3    2           2          8        8
# css_prob             10       NA        3   10          10          7       10

global_nagreements_goodmin_final <- median(global_nagreements_goodmin_opts)
id_nagreements_goodmin_final <- apply(id_nagreements_goodmin_opts, 2, median)


id_nagreements_goodmin_final[is.na(id_nagreements_goodmin_final)] <- global_nagreements_goodmin_final
id_nagreements_goodmin_final[id_nagreements_goodmin_final <= 2] <- global_nagreements_goodmin_final

print(id_nagreements_goodmin_final)
##
## Last filter: at least 5 agreements
##
id_nagreements_goodmin_final[id_nagreements_goodmin_final<5] <- 5
write.csv(id_nagreements_goodmin_final, file.path(dir_diag, "id_nagreements_goodmin_final.csv"))
write_rds(id_nagreements_goodmin_final, file.path(dir_diag, "id_nagreements_goodmin_final.rds"))

#------ Decide good/bad cells ------
for (id in ident_levels) {
    cli_h2(id)
}
df_css_list <- lapply(ident_levels, function(id) {
    cat(id, "... ")
    cutoff_id <- id_nagreements_goodmin_final[id]
    df_css_id <- df_css %>% dplyr::filter(css_id == id)
    is_good <- df_css_id$nagreements >= id_nagreements_goodmin_final[id]
    is_good[is.na(is_good)] <- FALSE

    css_id_refine <- as.character(df_css_id$css_id)
    css_id_refine[!is_good] <- "LOWCONF"
    df_css_id$css_id_refine <- css_id_refine

    return(df_css_id)
})
cat("\n")
df_css_refine <- do.call(rbind, df_css_list)
df_css_refine$css_id_refine <- standardize_factor(df_css_refine$css_id_refine, std_ident_levels)

df_css_refine <- df_css_refine[Cells(obj), ]
identical(Cells(obj), rownames(df_css_refine))
identical(Cells(obj), rownames(df_css))
write.csv(df_css_refine, file.path(dir_diag, "consensus.dataframe.csv"))
write_rds(df_css_refine, file.path(dir_diag, "consensus.dataframe.rds"))
write_parquet(df_css_refine, file.path(dir_diag, "consensus.dataframe.parquet"))

pal_ident2 <- c(pal_ident, "LOWCONF" = "lightgrey")

#------ barplot before and after QC ------
p1 <- ruok::qbarplot_table_cat(df_css_refine$css_id, do.prop.table = F, "css_id") + scale_fill_manual(values = pal_ident2)
p2 <- ruok::qbarplot_table_cat(df_css_refine$css_id_refine, do.prop.table = F, "css_id_refine") + scale_fill_manual(values = pal_ident2)
p <- wrap_plots(
    p1 + rremove("legend"),
    p2 + rremove("legend"),
    nrow = 1
)

ggsave(file.path(dir_diag, "barplot.pre_post_clean2.pdf"), p, width = 3, height = 4, useDingbats = F)
prop.table(table(df_css_refine$css_id_refine))

p1 <- qbarplot_table_catx(h = "css_id", v = "css_id_refine", df = df_css_refine, do.prop.table = F) + scale_fill_manual(values = pal_ident2)
p2 <- qbarplot_table_catx(h = "css_id", v = "css_id_refine", df = df_css_refine, do.prop.table = T) + scale_fill_manual(values = pal_ident2)
p <- wrap_plots(
    p1 + rremove("legend") + rotate_x_text(90) + labs(y = "num cells") + rremove("x.text"),
    p2 + rremove("legend") + rotate_x_text(90) + labs(y = "freq cells"),
    nrow = 2
)
ggsave(file.path(dir_diag, "barplot.clean2_each_ident.pdf"), p, width = 0.3 * length(ident_levels), height = 4.5, useDingbats = F)


#------ umap before and after QC ------
## use xenium.snippet_viz_categorical.R
dir_snippet_viz <- dir_diag
identical(Cells(obj), rownames(df_css_refine))
xmo <- AddMetaData(obj, df_css_refine)
pal_z <- pal_ident2
viz_what <- "css_id"
source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_categorical.R")

viz_what <- "css_id_refine"
source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_categorical.R")


#------ dotplot before and after QC ------
## use std.dotplot_pretty.cellstates.R

#------------------ ~~~ Export object to clean2 ~~~ --------------------
cli_h1("Export object to clean2")
cat("subsetting...")
xmo <- subset(xmo, css_id_refine != "LOWCONF")
cat("done\n")
stopifnot(identical(as.character(xmo$css_id), as.character(xmo$css_id_refine)))
xmo$css_id_refine <- NULL

cat("exporting...")
write_seurat(xmo, dir_obj, "pass")
cat("done\n")

#------ color legend ------
ggsave(file.path(dir_diag, "legend.css_id_refine.pdf"),
    pal_to_ggplot(pal_ident2, "css_id_refine"),
    width = 3, height = 5
)
ggsave(file.path(dir_diag, "legend.css_id.pdf"),
    pal_to_ggplot(pal_ident, "css_id"),
    width = 3, height = 5
)

cat("[done]")
timestamp()
