suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Consensus decision of MapQuery's predicted.id (`cat_transfer`)
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: Nov 10, 2024
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
    library(ruok)
    library(SCpubr)
    library(ComplexHeatmap)
    library(patchwork)
    my_scatter_themevoid <- theme_pubr(base_size = 6, legend = "right") %+replace% theme(
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = rel(1)),
        axis.line = element_blank()
    )

    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.nmf.viz.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
})
cmdargs <- commandArgs(trailingOnly = TRUE)

if (length(cmdargs) > 0) {
    f_sp <- cmdargs[1]
    cat_transfer <- cmdargs[2] # celltype | cell_state_paper
} else {
    # f_sp <- "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_xenium5k/celltype_T/integrate_seurat_50pc/ready.sr3.rds"
    # cat_transfer <- "cell_state_paper"

    f_sp <- "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_xenium5kPlus/celltype_B/integrate_seurat_50pc/ready.sr3.rds"
    cat_transfer <- "cell_state_paper"
    f_sp <- "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_xenium5kPlus/celltype_Endo/integrate_seurat_50pc/ready.sr3.rds"
    cat_transfer <- "cell_state_paper"
}

#------------------ ~~~ Env variables/helpers ~~~ --------------------
cli_h1("Env variables/helpers")
sc_assay <- "XNA"
sp_assay <- "integrated"
mapquery_method <- "cca"
anchor_method <- transfer_method <- mapquery_method

n_gene_portions <- 10 ## do not change
#------------------ ~~~ Setup outdirs ~~~ --------------------
dir_res <- file.path(
    dirname(f_sp),
    sprintf("mapquery_css_%s.ref_%s.query_%s", anchor_method, sc_assay, sp_assay),
    cat_transfer
)
fs::dir_create(dir_res)

#------------------ ~~~ Read in  ~~~ --------------------
cli_h1("Read-in")
#------ query object ------
sp <- read_rds(f_sp)
DefaultAssay(sp) <- sp_assay
print(sp)

#------ results of the 10 gene portions ------
dir_anchor_portions_prefix <- file.path(
    dirname(f_sp),
    sprintf(
        "mapquery_worker_%s.ref_%s.query_%s.gene_portion_%s",
        anchor_method, sc_assay, sp_assay, 1:n_gene_portions
    )
)
names(dir_anchor_portions_prefix) <- 1:n_gene_portions

df_pred_list <- lapply(names(dir_anchor_portions_prefix), function(x) {
    df <- read_rds(file.path(
        dir_anchor_portions_prefix[x],
        mapquery_method, cat_transfer, "df.TransferData.rds"
    ))
    df <- rownames_to_column(df, "cellname")
    df$portion <- x
    df$sample <- str_extract(df$cellname, "ART[0-9]+")
    return(df)
})

mat_pred_list <- lapply(names(dir_anchor_portions_prefix), function(x) {
    mat <- read_rds(file.path(
        dir_anchor_portions_prefix[x],
        mapquery_method, cat_transfer, "pred_prob.mat.rds"
    ))
    mat
})

#------ load std levels ------
std_ident_levels <- read_rds(
    file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate",
        "pat102/atlas",
        "idents_levels_cellstates.rds"
    )
)
#------ setup colors ------
cnames <- Cells(sp)
ident_levels <- Reduce(unique, lapply(mat_pred_list, rownames))
ident_levels <- intersect(std_ident_levels, ident_levels)
print(ident_levels)
pal_ident <- init_pal_d(ident_levels)
pal_ident2 <- c(pal_ident, `LOWCONF` = "lightgrey")
ggsave(file.path(dir_res, sprintf("legend.%s.pdf", cat_transfer)),
    pal_to_ggplot(pal_ident, cat_transfer),
    width = 4, height = 4, useDingbats = F
)
ggsave(file.path(dir_res, sprintf("legend.%s.2.pdf", cat_transfer)),
    pal_to_ggplot(pal_ident2, cat_transfer),
    width = 4, height = 4, useDingbats = F
)

#------------------ ~~~ Consensus identities ~~~ --------------------
cli_h1("Consensus identities")
## Assess if the predicted.id are consisitent across portions
## a long df: cell | portion_id | predicted.id
df_pred <- do.call("rbind", df_pred_list)
df_pred$predicted.id <- standardize_factor(df_pred$predicted.id, std_ident_levels)
df_pred$portion <- as.factor(as.numeric(df_pred$portion))
tail(df_pred)
#                 cellname predicted.id.score predicted.id portion sample
# 471795 ART304_ojakciea-1          0.3513797      CD4-TCM      10 ART304
# 471796 ART304_ojancgeo-1          0.3457251     CD4-TREG      10 ART304
# 471797 ART304_ojanfhoc-1          0.4617882      CD8-TRM      10 ART304
pal_sample <- init_pal_d(df_pred$sample, pal = "parade")

p <- df_pred %>%
    ggplot(., aes(x = portion, fill = predicted.id)) +
    geom_bar() +
    scale_fill_manual(values = pal_ident) +
    labs(y = "num cells")
p <- p + labs(caption = report_chisq(
    chisq.test(table(df_pred$predicted.id, df_pred$portion))
))

ggsave(file.path(dir_res, "barplot.ident_by_portions.pdf"),
    p + rremove("legend"),
    width = 4, height = 3, useDingbats = F
)

#------ Determine the consensus idents ------
f_df_css <- file.path(dir_res, "consensus.dataframe.rds")
df_css_opts <- df_pred %>%
    dplyr::group_by(cellname, predicted.id) %>%
    dplyr::summarise(
        nportion = n(),
        predicted.id.score = mean(predicted.id.score, na.rm = TRUE)
    )

# if (!file.exists(f_df_css)) {
if (TRUE) {
    df_css_res <- df_css_opts %>%
        dplyr::group_by(cellname) %>%
        dplyr::arrange(desc(nportion), desc(predicted.id.score), .by_group = TRUE) %>%
        dplyr::slice_head(n = 1)
    print(head(df_css_res)) # ! export
    # cellname         predicted.id nportion predicted.id.score
    #   <chr>            <chr>           <int>              <dbl>
    # 1 ART10_aabghhkd-1 CD4-TCM             2              0.373
    # 2 ART10_aabgophb-1 CD4-TN              9              0.478
    # 3 ART10_aabigcah-1 CD4-TCM             7              0.467
    colnames(df_css_res) <- c("cellname", "css_id", "nagreements", "css_prob")

    if (!identical(df_css_res$cellname, Cells(sp))) {
        message("unifying cells...")
        idx <- match(Cells(sp), df_css_res$cellname)
        df_css_res <- df_css_res[idx, ]
        identical(df_css_res$cellname, Cells(sp))
        df_css_res$cellname <- Cells(sp)
        df_css_res$css_id <- standardize_factor(df_css_res$css_id, std_ident_levels)
    }
    df_css_res <- as.data.frame(df_css_res)
    rownames(df_css_res) <- df_css_res$cellname
    
    ## export
    write_csv(df_css_res, file.path(dir_res, "consensus.dataframe.csv"))
    write_rds(df_css_res, file.path(dir_res, "consensus.dataframe.rds"))
    write_parquet(df_css_res, file.path(dir_res, "consensus.dataframe.parquet"))
} else {
    df_css_res <- read_rds(f_df_css)
}

#------------------ ~~~ Consensus probs ~~~ --------------------
cli_h1("Consensus probs")
f_mat_prob <- file.path(dir_res, "pred_prob_mean.mat.rds")
# if (!file.exists(f_mat_prob)) {
if (TRUE) {

    mat_pred_list <- lapply(mat_pred_list, function(mat) {
        mat[ident_levels, cnames]
    })
    array_pred <- simplify2array(mat_pred_list)
    mat_prob <- apply(array_pred, c(1, 2), mean, na.rm = TRUE)
    if (!identical(colnames(mat_prob), Cells(sp))) {
        message("unifying cells...")
        idx <- match(Cells(sp), colnames(mat_prob))
        mat_prob <- mat_prob[, idx]
    }
    mat_prob <- mat_prob[ident_levels, ]
    write_rds(mat_prob, f_mat_prob)
} else {
    mat_prob <- read_rds(f_mat_prob)
}


#------------------ ~~~ Diagnosis consensus identities ~~~ --------------------
cli_h1("Diagnoise css_id")
head(df_css_res)
stopifnot(identical(Cells(sp), rownames(df_css_res)))
stopifnot(identical(colnames(mat_prob), rownames(df_css_res)))

#------ heatmap of prob or nAgreements ------
head(df_css_res)
head(df_css_opts)
set.seed(1026)
n_demo_cells_each <- 100
demo_cnames <- df_css_res %>%
    dplyr::group_by(css_id) %>%
    dplyr::slice_sample(n = n_demo_cells_each) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(css_id, desc(nagreements), desc(css_prob))
head(demo_cnames)
demo_cnames <- demo_cnames %>% dplyr::pull(cellname)
str(demo_cnames)

demo_mat_css <- df_css_opts %>%
    dplyr::filter(cellname %in% demo_cnames) %>%
    tidyr::pivot_wider(., id_cols = cellname, names_from = predicted.id, values_from = nportion, values_fill = 0)
demo_mat_css <- column_to_rownames(demo_mat_css, "cellname")
demo_mat_css <- as.matrix(demo_mat_css)
demo_mat_css <- demo_mat_css[, intersect(ident_levels, colnames(demo_mat_css))]
demo_mat_css <- demo_mat_css[demo_cnames, ]

demo_mat_prob <- t(mat_prob[, demo_cnames])

demo_cnames_to_ident <- df_css_res %>%
    dplyr::filter(cellname %in% demo_cnames) %>%
    dplyr::select(cellname, css_id) %>%
    deframe()
demo_cnames_to_ident <- demo_cnames_to_ident[demo_cnames]

demo_cnames_htmp_anno <- rowAnnotation(
    cluster = demo_cnames_to_ident, col = list(cluster = pal_ident),
    show_annotation_name = c(FALSE)
)
demo_ident_htmp_anno <- columnAnnotation(
    cluster = ident_levels,
    col = list(cluster = pal_ident),
    show_legend = c(FALSE), show_annotation_name = c(FALSE)
)

adhoc_heatmap_opt_list <- list(
    show_row_names = F, cluster_rows = F,
    row_title = sprintf("%d cells per cluster", n_demo_cells_each),
    cluster_columns = F, show_column_names = F, column_title = "prediction",
    left_annotation = demo_cnames_htmp_anno,
    top_annotation = demo_ident_htmp_anno,
    use_raster = T, raster_by_magick = TRUE
)

pdf(file.path(dir_res, "heatmap.demo.nAgreements.pdf"), width = 3.5, height = 5, useDingbats = F)
p <- do.call(Heatmap, c(
    list(
        matrix = demo_mat_css, name = "nagreements",
        col = heatmap_color_fun_cont(1, 9, "Purples 2", rev = T)
    ),
    adhoc_heatmap_opt_list
))
draw(p)
dev.off()
pdf(file.path(dir_res, "heatmap.demo.pred_prob.pdf"), width = 3.5, height = 5, useDingbats = F)
p <- do.call(Heatmap, c(
    list(
        matrix = demo_mat_prob, name = "prob",
        col = heatmap_color_fun_cont(0.1, 0.9, "Purples 2", rev = T)
    ),
    adhoc_heatmap_opt_list
))
draw(p)
dev.off()

#------ hist of prob or nAgrerments ------
head(df_css_res)
library(ggridges)
for (tmp in c("nagreements", "css_prob")) {
    p <- ggplot(df_css_res, aes_string(x = tmp, y = "css_id", fill = "css_id")) +
        geom_density_ridges_gradient(scale = 0.6, rel_min_height = 0.01) +
        scale_fill_manual(values = pal_ident) +
        scale_y_discrete(limits = rev)
    ggsave(file.path(dir_res, sprintf("histogram.%s.pdf", tmp)),
        p + rremove("legend"),
        width = 4, height = 0.5 * length(ident_levels), useDingbats = FALSE
    )
}

#------ prob/nUMI/nCount vs nAgreements for all cells ------
stopifnot(identical(rownames(df_css_res), Cells(sp)))
df_css_res$nCount_Xenium <- sp$nCount_Xenium
df_css_res$nFeature_Xenium <- sp$nFeature_Xenium
df_css_res$sample <- sp$sample
df_css_res$sample <- standardize_factor(df_css_res$sample, names(pal_sample))

for (basic_metric in c("nCount_Xenium", "nFeature_Xenium", "css_prob")) {
    cat(basic_metric, "... ")
    p <- df_css_res %>%
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
            dir_res,
            sprintf("violin.%s.all_cells.pdf", basic_metric)
        ),
        p,
        width = 4, height = 2.5, useDingbats = F
    )
}
cat("\n")

#------ prob/nUMI/nCount vs nAgreements faceted by cell ident  ------
head(df_css_res)
for (basic_metric in c("nCount_Xenium", "nFeature_Xenium", "css_prob")) {
    cat(basic_metric, ">>> ")
    fs::dir_create(
        file.path(dir_res, sprintf("violin.%s.faceted_by_ident", basic_metric))
    )

    for (id in ident_levels) {
        cat(id, "...")
        p <- df_css_res %>%
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
                file.path(dir_res, sprintf("violin.%s.faceted_by_ident", basic_metric)),
                sprintf("violin.%s.pdf", id)
            ),
            p,
            width = 4, height = 2.5, useDingbats = F
        )
    }
    cat("\n")
}
cat("\n")

#------ export ------
write_csv(df_css_res, file.path(dir_res, "consensus.dataframe.csv"))
write_rds(df_css_res, file.path(dir_res, "consensus.dataframe.rds"))
write_parquet(df_css_res, file.path(dir_res, "consensus.dataframe.parquet"))

#------------------ ~~~ Consensus UMAPs ~~~ --------------------
f_css_umap <- file.path(dirname(dir_res), "ref.umap.projection.dataframe.rds")
# if (file.exists(f_css_umap)) {
if (FALSE){
    umap_css <- read_rds(f_css_umap)
} else {
    umap_list <- lapply(names(dir_anchor_portions_prefix), function(x) {
        umap <- read_rds(file.path(
            dir_anchor_portions_prefix[x],
            "ref.umap.projection.dataframe.rds"
        ))
        umap <- umap[Cells(sp), ]
        return(umap)
    })
    array_umap <- simplify2array(umap_list)
    umap_css <- apply(array_umap, c(1, 2), median, na.rm = TRUE)
    stopifnot(identical(rownames(umap_css), Cells(sp)))
    write_rds(umap_css, f_css_umap)
}
umap_css <- CreateDimReducObject(
    embeddings = umap_css,
    assay = DefaultAssay(sp)
)
sp[["ref.umap"]] <- umap_css
print(sp)

#------------------ ~~~ Visualization ~~~ --------------------
cli_h1("Visualization")
sp <- AddMetaData(
    sp,
    metadata = df_css_res[, setdiff(colnames(df_css_res), colnames(sp@meta.data))]
)

## single umap showing css_id
for (umap_use in c("umap", "ref.umap")) {
    p <- DimPlot(sp,
        group.by = "css_id",
        reduction = umap_use,
        pt.size = 3, raster = T, raster.dpi = c(1024, 1024),
        shuffle = T, label = TRUE, cols = pal_ident2
    ) + my_scatter_themevoid
    ggsave(file.path(dir_res, sprintf("dimplot.%s.%s.pdf", umap_use, "css_id")),
        p + rremove("legend"),
        width = 6, height = 6, useDingbats = F
    )
}
## umap facted by css_id
for (umap_use in c("umap", "ref.umap")) {
    cat(umap_use, ">>> ")
    dir_create(file.path(dir_res, sprintf("dimplot.%s.%s.faceted", umap_use, "css_id")))
    pdf(
        file.path(
            dir_res, sprintf("dimplot.%s.%s.faceted", umap_use, "css_id"),
            "dimplot_facted.%02d.pdf"
        ),
        width = 6, height = 6, useDingbats = F, onefile = FALSE
    )
    p_list <- list()
    for (id in ident_levels) {
        cat(id, "... ")

        if (! id %in% sp@meta.data[["css_id"]]) {
            cat("no cells present...")
            next()
        }

        p <- DimPlot(sp,
            cells.highlight = Cells(sp)[sp[["css_id"]] == id],
            reduction = umap_use,
            pt.size = 3, raster = T, raster.dpi = c(1024, 1024),
            sizes.highlight = 3,
            label = FALSE, cols.highlight = unname(pal_ident2[id]),
        ) + my_scatter_themevoid
        p <- p + rremove("legend") + labs(title = id)
        print(p)
        p_list <- c(p_list, list(p))
    }
    dev.off()

    p_combo <- wrap_plots(p_list, ncol = ceiling(sqrt(length(p_list))))
    ggsave(
        file.path(
            dir_res, sprintf("dimplot.%s.%s.faceted_combo.pdf", umap_use, "css_id")
        ),
        p_combo,
        limitsize = FALSE,
        width = 3 * ceiling(sqrt(length(p_list))),
        height = 3 * ceiling(sqrt(length(p_list))), useDingbats = F
    )

    cat("\n")
}
#------------------ ~~~ Export to XeniumExplorer format ~~~ --------------------

for (ident_str in c("css_id")) {
    df_css_res %>%
        df_to_exnium_explorer(., "cellname", ident_str) %>%
        write_csv(., file.path(
            dir_res,
            sprintf("to_xenium_explorer_groups.%s.csv", ident_str)
        ))
}



cli_alert_success("[done]")
timestamp()
