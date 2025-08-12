suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: CCA-based MapQuery
# This is a worker in part of the Consensus MapQuery workflow.
# It performs MapQuery using 9 portions of genes and leaves the 1 portion for evaluation.
#
# It also TransferData of the expression of the 1 evalution gene portion.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: Nov 8, 2024
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
    library(scattermore)
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
})
cmdargs <- commandArgs(trailingOnly = TRUE)
print(cmdargs)

if (length(cmdargs) > 0) {
    f_sc <- cmdargs[[1]]
    f_sp <- cmdargs[[2]]
    cat_transfer <- cmdargs[[3]] # celltype | cell_state_paper
    sc_assay <- cmdargs[[4]] # RNA|XNA
    sp_assay <- cmdargs[[5]] # Xenium|integrated
    dir_gene_portions <- cmdargs[[6]] # xenium5k | xenium5kPlus
    gene_portion_index <- cmdargs[[7]] # 1,2,...10
    gene_portion_index <- as.numeric(gene_portion_index)
} else {
    # f_sc <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/T/downto_1000_by_cell_state_paper/test_power_xenium5k/ready.seurat.rds"
    # f_sp <- "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_xenium5k/celltype_T/integrate_seurat_50pc/ready.sr3.rds"

    # cat_transfer <- "cell_state_paper"

    # sc_assay <- "XNA"
    # sp_assay <- "integrated"

    # gene_portion_index <- 4
    # dir_gene_portions <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/T/downto_1000_by_cell_state_paper/gene_portions_xenium5k"
}
#------------------ ~~~ Env helpers/parameters ~~~ --------------------
param.FindTransferAnchors_n_dims <- 50 # old: 30
mapquery_method <- "cca"
anchor_method <- "cca"
transfer_method <- mapquery_method

ref_umap_name <- "umap.xna"

do_force_findanchor <- TRUE
do_force_transfergene <-  TRUE
do_force_mapquery <- TRUE

#------------------ ~~~ Output ~~~ --------------------
dir_anchor <- file.path(
    dirname(f_sp),
    sprintf(
        "mapquery_worker_%s.ref_%s.query_%s.gene_portion_%s",
        anchor_method, sc_assay, sp_assay, gene_portion_index
    )
)
dir_cat_transfer <- file.path(dir_anchor, transfer_method, cat_transfer)
dir_anchor_eval <- file.path(dir_anchor, transfer_method, "eval_anchors")
fs::dir_create(dir_anchor)
fs::dir_create(dir_cat_transfer)
fs::dir_create(dir_anchor_eval)


#------------------ ~~~ Input ~~~ --------------------
sc <- read_rds(f_sc)
sp <- read_rds(f_sp)

DefaultAssay(sc) <- sc_assay
DefaultAssay(sp) <- sp_assay
print(sc)
print(sp)
adhoc_load_gene_portions <- function(
    dir_gene_portions, eval_gene_portion_index, n_portions = 10) {
    f_test <- file.path(
        dir_gene_portions,
        sprintf("genes_of_portion_%d.rds", eval_gene_portion_index)
    )
    f_train_list <- file.path(
        dir_gene_portions,
        sprintf(
            "genes_of_portion_%d.rds",
            setdiff(1:n_portions, eval_gene_portion_index)
        )
    )
    g_test <- read_rds(f_test)
    g_train <- as.character(unlist(lapply(f_train_list, read_rds)))
    return(list(
        genes_test = g_test,
        genes_train = g_train
    ))
}

## read genes for training and testing
genes_combo <- adhoc_load_gene_portions(
    dir_gene_portions = dir_gene_portions,
    eval_gene_portion_index = gene_portion_index
)
str(genes_combo)


#------------------ ~~~ FindAnchors ~~~ --------------------
cli_h1("FindTransferAnchors")
timestamp()
f_anchors <- file.path(dir_anchor, "anchors.rds")

## Make sure the sc ref use the XNA related reductions
sc[["umap"]] <- NULL
sc[["pca"]] <- NULL
sc[["pca"]] <- sc[["pca.xna"]]
sc[["umap"]] <- sc[["umap.xna"]]

if ( (!file.exists(f_anchors)) | do_force_findanchor ) {
    anchors <- FindTransferAnchors(
        reference = sc, query = sp,
        reference.assay = sc_assay, query.assay = sp_assay,
        features = genes_combo$genes_train,
        reduction = "cca",
        k.filter = NA,
        dims = 1:param.FindTransferAnchors_n_dims
    )
    write_rds(anchors, f_anchors)
    cli_alert_success("[finished] FindTransferAnchors.")
} else {
    cli_alert_success("reading the existing anchors...")
    anchors <- read_rds(f_anchors)
}
timestamp()

#------------------ ~~~ TransferData genes ~~~ --------------------
cli_h1("TransferData - genes")
timestamp()
f_impute_data <- file.path(dir_anchor_eval, "TransferData_impute_genes.matrix.rds")
if ( (!file.exists(f_impute_data)) | do_force_transfergene ) {
    refdata <- GetAssayData(sc, slot = "data", assay = sc_assay)
    imputation <- TransferData(
        anchorset = anchors,
        query = sp,
        refdata = refdata,
        weight.reduction = "cca",
        dims = 1:param.FindTransferAnchors_n_dims,
        slot = "data"
    )
    class(imputation) # seurat
    print(imputation@assays$id) # the new/imputed assay gene data
    print(dim(imputation@assays$id@data))

    impute_genes_all <- imputation@assays$id@data
    write_rds(imputation@assays$id@data, f_impute_data)
    cli_alert_success("[finished] TransferData")
} else {
    cli_alert_success("reading the existing TransferData...")
    impute_genes_all <- read_rds(f_impute_data)
}
timestamp()

#------ Fetch the imputed & observed gene expressions ------
## to prepare for future evalutions
## The observed gene expressions fetch the assay (i.e. `sp_assay`) used for finding anchor.
## If the sp's assay is 'integrated', xenium is not used for finding anchors,
## so it is unfair to fetch the observed expression from xenium assay.

f_sp_impute_g_tain <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_imputed_genes_train"))
f_sp_impute_g_test <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_imputed_genes_test"))
f_sp_observe_g_train <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_observed_genes_train"))
f_sp_observe_g_test <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_observed_genes_test"))

if ( !file_exists(f_sp_observe_g_test) | do_force_transfergene ) {
    ## impute_sp_g_train + obs_sp_g_train: evaluating training performance
    ## impute_sp_g_test  + obs_sp_g_test : evaluating testing performance

    obs_genes_all <- GetAssayData(sp, slot = "data", assay = sp_assay)

    genes_shared <- intersect(rownames(impute_genes_all), rownames(obs_genes_all))
    str(genes_shared)
    genes_combo <- lapply(genes_combo, function(x) intersect(x, genes_shared))

    impute_sp_g_train <- impute_genes_all[genes_combo$genes_train, ]
    impute_sp_g_test <- impute_genes_all[genes_combo$genes_test, ]
    obs_sp_g_train <- obs_genes_all[genes_combo$genes_train, ]
    obs_sp_g_test <- obs_genes_all[genes_combo$genes_test, ]

    stopifnot(identical(colnames(impute_sp_g_train), Cells(sp)))
    stopifnot(identical(colnames(impute_sp_g_test), Cells(sp)))
    stopifnot(identical(colnames(obs_sp_g_train), Cells(sp)))
    stopifnot(identical(colnames(obs_sp_g_test), Cells(sp)))

    write_rds(impute_sp_g_train, f_sp_impute_g_tain)
    write_rds(impute_sp_g_test, f_sp_impute_g_test)
    write_rds(obs_sp_g_train, f_sp_observe_g_train)
    write_rds(obs_sp_g_test, f_sp_observe_g_test)
} else {
    impute_sp_g_train <- read_rds(f_sp_impute_g_tain)
    impute_sp_g_test <- read_rds(f_sp_impute_g_test)
    obs_sp_g_train <- read_rds(f_sp_observe_g_train)
    obs_sp_g_test <- read_rds(f_sp_observe_g_test)
}


#------------------ ~~~ MapQuery ~~~ --------------------
cli_h1("MapQuery")
timestamp()
f_TransferData <- file.path(
    dir_cat_transfer,
    sprintf("df.TransferData.rds")
)
xmo_old_umap <- sp[["umap"]]
meta_colnames_old <- colnames(sp@meta.data)
umap_mapquery_name <- "ref.umap"
if ( (!file.exists(f_TransferData)) | do_force_mapquery ) {
    sp <- MapQuery(
        anchorset = anchors,
        reference = sc, query = sp,
        refdata = sc@meta.data[[cat_transfer]],
        reduction.model = ref_umap_name,
        transferdata.args = list(slot = "data")
    )
    cat("These are the new added data:")
    print(setdiff(colnames(sp@meta.data), meta_colnames_old))
    print(sp)
    ## a new umap called 'ref.umap'
    ## a new reduction 'cca' is created.
    ## the imputed data is missing
    ##

    umap_mapquery_name <- "ref.umap"
    pred <- sp@meta.data[, setdiff(colnames(sp@meta.data), meta_colnames_old)]
    head(pred)

    ## export predicted labels
    write_rds(pred, f_TransferData)
    write.csv(pred, paste0(f_TransferData, ".csv"))
    nanoparquet::write_parquet(pred, paste0(f_TransferData, ".parquet"))

    ## export the UMAP space of sc such that sp --projected--> sc
    write_rds(
        sp[[umap_mapquery_name]],
        file.path(dir_anchor, "ref.umap.projection.seurat.rds")
    )
    write_rds(
        Embeddings(sp, umap_mapquery_name),
        file.path(dir_anchor, "ref.umap.projection.dataframe.rds")
    )

    ## export the prediction prob for each cell identity
    mat_pred_prob <- GetAssayData(sp, slot = "data", assay = "prediction.score.id")
    write_rds(
        mat_pred_prob,
        file.path(dir_cat_transfer, "pred_prob.mat.rds")
    )

    ## export the updated object
    try(write_seurat(sp, dir_anchor, "ready"))
} else {
    pred <- read_rds(f_TransferData)
    print(head(pred))
    sp <- read_rds(file.path(dir_anchor, "ready.seurat.rds"))
    mat_pred_prob <- read_rds(file.path(dir_cat_transfer, "pred_prob.mat.rds"))
}

#------------------ ~~~ Basic visualization ~~~ --------------------
std_ident_levels <- read_rds(
    file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate",
        "pat102/atlas",
        "idents_levels_cellstates.rds"
    )
)

sp$predicted.id <- standardize_factor(sp$predicted.id, std_ident_levels)
ident_levels <- levels(sp$predicted.id)
print(ident_levels)
Idents(sp) <- "predicted.id"
print(table(Idents(sp)))
pred$predicted.id <- standardize_factor(pred$predicted.id, std_ident_levels)

cell_type_compartment_name <- str_extract(f_sp, "celltype_[A-Za-z]+") %>% str_remove_all("celltype_")
# if (cell_type_compartment_name %in% c("T", "B", "Mye")) {
#     pal_use <- init_pal_d(pred$predicted.id)
# } else {
#     pal_use <- init_pal_d(pred$predicted.id, pal = "circus")
# }
## use the same color pal for cell states
pal_use <- init_pal_d(pred$predicted.id)

pl <- enframe(pal_use) %>%
    mutate(name = factor(name, levels = name)) %>%
    ggplot() +
    geom_point(aes(x = name, y = factor(1), color = name), size = 3) +
    scale_color_manual(values = pal_use) +
    labs(color = cat_transfer)
ggsave(file.path(dir_cat_transfer, sprintf("legend.%s.pdf", cat_transfer)),
    as_ggplot(get_legend(pl)),
    width = 3, height = 3, useDingbats = F
)
rm(pl)

#------ prob hist ------

gghistogram(pred, x = "predicted.id.score", bins = 100)
library(ggridges)
p <- ggplot(pred, aes(x = predicted.id.score, y = predicted.id, fill = predicted.id)) +
    geom_density_ridges_gradient(scale = 0.6, rel_min_height = 0.01) +
    scale_fill_manual(values = pal_use) +
    scale_y_discrete(limits = rev)
ggsave(file.path(dir_cat_transfer, "hist_prob.pdf"),
    p + rremove("legend"),
    width = 4, height = 0.5 * length(ident_levels), useDingbats = FALSE
)

#------ heatmap prob ------
dim(mat_pred_prob)
ident_levels_with_prob <- intersect(ident_levels, rownames(mat_pred_prob))
mat_pred_prob <- mat_pred_prob[ident_levels_with_prob, , drop = FALSE]

set.seed(1026)
pred_dns <- pred %>%
    tibble::rownames_to_column("cname") %>%
    dplyr::group_by(predicted.id) %>%
    dplyr::slice_sample(n = 100) %>%
    dplyr::arrange(desc(predicted.id.score), .by_group = TRUE)
table(pred_dns$predicted.id)
head(pred_dns)
tmp <- t(as.matrix(mat_pred_prob[, pred_dns$cname]))

p <- Heatmap(
    tmp,
    col = heatmap_color_fun_cont(0.1, 0.9, "Purples 2", rev=T),
    name = "prob",
    show_row_names = F, cluster_rows = F, row_title = "cells",
    cluster_columns = F, show_column_names = F, column_title = "prediction",
    left_annotation = rowAnnotation(
        cluster = pred_dns$predicted.id,
        col = list(cluster = pal_use), 
        show_annotation_name = c(FALSE)
    ),
    top_annotation = columnAnnotation(
        cluster = colnames(tmp), 
        col =  list(cluster = pal_use), 
        show_legend = c(FALSE), 
        show_annotation_name = c(FALSE)
    ),
    use_raster = T, raster_by_magick = TRUE
)
pdf(
    file.path(
        dir_cat_transfer,
        sprintf("heatmap_prob.%s.pdf", cat_transfer)
    ),
    width = 3.5, height = 5, useDingbats = F
)
draw(p)
dev.off()



#------ barplot ------
p <- qbarplot_table_cat(pred$predicted.id, name_x = cat_transfer, do.prop.table = FALSE) +
    scale_fill_manual(values = pal_use)
ggsave(file.path(dir_cat_transfer, "barplot.count.pdf"),
    p + rremove("legend"),
    width = 2, height = 3, useDingbats = FALSE
)
p <- qbarplot_table_cat(pred$predicted.id, name_x = cat_transfer, do.prop.table = TRUE) +
    scale_fill_manual(values = pal_use)
ggsave(file.path(dir_cat_transfer, "barplot.freq.pdf"),
    p + rremove("legend"),
    width = 2, height = 3, useDingbats = FALSE
)
#------ umap ------
for (umap_use in c("umap", umap_mapquery_name)) {
    p <- DimPlot(sp,
        group.by = "predicted.id",
        reduction = umap_use,
        pt.size = 3, raster = T, raster.dpi = c(1024, 1024),
        shuffle = T, label = TRUE, cols = pal_use
    ) + my_scatter_themevoid
    ggsave(file.path(dir_cat_transfer, sprintf("dimplot.%s.%s.pdf", umap_use, "predicted.id")),
        p + rremove("legend"),
        width = 6, height = 6, useDingbats = F
    )
}
for (umap_use in c("umap", umap_mapquery_name)) {
    cat(umap_use, ">>> ")
    dir_create(file.path(dir_cat_transfer, sprintf("dimplot.%s.%s.faceted", umap_use, "predicted.id")))
    pdf(
        file.path(
            dir_cat_transfer, sprintf("dimplot.%s.%s.faceted", umap_use, "predicted.id"),
            "dimplot_facted.%02d.pdf"
        ),
        width = 6, height = 6, useDingbats = F, onefile = FALSE
    )
    p_list <- list()
    for (id in ident_levels) {
        cat(id, "... ")
        p <- DimPlot(sp,
            cells.highlight = Cells(sp)[Idents(sp) == id],
            reduction = umap_use,
            pt.size = 3, raster = T, raster.dpi = c(1024, 1024),
            sizes.highlight = 3,
            label = FALSE, cols.highlight = unname(pal_use[id]),
        ) + my_scatter_themevoid
        p <- p + rremove("legend") + labs(title = id)
        print(p)
        p_list <- c(p_list, list(p))
    }
    dev.off()

    p_combo <- wrap_plots(p_list, ncol = ceiling(sqrt(length(p_list))))
    ggsave(
        file.path(
            dir_cat_transfer, sprintf("dimplot.%s.%s.faceted_combo.pdf", umap_use, "predicted.id")
        ),
        p_combo,
        limitsize = FALSE,
        width = 3 * ceiling(sqrt(length(p_list))),
        height = 3 * ceiling(sqrt(length(p_list))), useDingbats = F
    )

    cat("\n")
}


cat("[done]")
timestamp()
