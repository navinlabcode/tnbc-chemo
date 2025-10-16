suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Identify niches using the cell composition data frames of ROIs
#
# Input: data frame of cell type compositions of ROIs. I work with either all samples or one sample at a time.
#
# Ref: https://github.com/satijalab/seurat/blob/e44cb2ce21aff3cbd94c2ddf24f9a1db6745e121/R/utilities.R#L3036
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025 Aug 26
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
    library(future)
    options(future.globals.maxSize = 8 * 1024^3) # 8GB
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    dir_compo <- cmdargs[1]
    f_cell_meta <- cmdargs[2]
} else {
    obj_suit <- "ALL"
    radius_final_use <- 30
    f_cell_meta <- file.path(
        "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/",
        obj_suit, "inputs", 
        "dataframe.cellmeta_combo.rds"
    )
    dir_compo <- file.path(
        "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/",
        obj_suit, "scNicheRadius", glue("R{radius_final_use}")
    )
    sample <- 'ART276'
    dir_compo=glue("/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/ALL/scNicheRadius/R30/split_sample/{sample}")
    f_cell_meta=glue("/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/ALL/inputs/split_sample/{sample}/inputs/dataframe.cellmeta_combo.rds")
}

#------ Outputs ------
dir_res <- file.path(dir_compo, "niche_proposal_onestep")
fs::dir_create(dir_res)

#------------------ ~~~ Load Inputs ~~~ --------------------
cli_h1("Load Inputs")

f_niche_composition_freq <- file.path(dir_compo, "result.full.composition_roi_freq.mat.rds")
f_niche_composition_ncells <- file.path(dir_compo, "result.full.composition_roi_ncell.mat.rds")

mat_niche_compo_freq <- read_rds(f_niche_composition_freq)
mat_niche_compo_ncells <- read_rds(f_niche_composition_ncells)
df_cell_meta <- read_rds(f_cell_meta)
print(head(df_cell_meta))
print(head(mat_niche_compo_freq))
# stopifnot(sum(mat_niche_compo_freq[1, ]) == 1)
if (!identical(rownames(mat_niche_compo_freq), rownames(df_cell_meta))) {
    shared_cells <- intersect(rownames(mat_niche_compo_freq), rownames(df_cell_meta))
    mat_niche_compo_freq <- mat_niche_compo_freq[shared_cells, , drop = FALSE]
    mat_niche_compo_ncells <- mat_niche_compo_ncells[shared_cells, , drop = FALSE]
    df_cell_meta <- df_cell_meta[shared_cells, , drop = FALSE]
    str(shared_cells)
}
identical(colnames(mat_niche_compo_freq), colnames(mat_niche_compo_ncells))
feature_names <- colnames(mat_niche_compo_freq)
feature_names <- str_replace_all(feature_names, "module_score_", "")
str(feature_names)
colnames(mat_niche_compo_freq) <- colnames(mat_niche_compo_ncells) <- feature_names

#------------------ ~~~ Initiate Niche Assay ~~~ --------------------
cli_h1("Initiate Niche Assay")
f_niche_obj <- file.path(dir_res, "niche_sr_object.rds")
# if (file.exists(f_niche_obj)) {
if (F) {
    niche_assay <- read_rds(f_niche_obj)
    print(niche_assay)
} else {
    niche_assay <- CreateSeuratObject(
        counts = t(mat_niche_compo_ncells),
        assay = "niche",
        meta.data = df_cell_meta[rownames(mat_niche_compo_ncells), , drop = FALSE], 
        min.cells = 1, min.features = 1)
    print(niche_assay)
    write_rds(niche_assay, file.path(dir_res, "niche_sr_object.rds"))

    niche_assay <- NormalizeData(niche_assay, assay = "niche", 
        normalization.method = "RC", 
        scale.factor = 100)
    write_rds(niche_assay, file.path(dir_res, "niche_sr_object.rds"))
}

min_ncount_niche <- min(niche_assay@meta.data[, 'nCount_niche'])
if (min_ncount_niche < 3) {
    cli_alert_danger("Some ROIs have very low nCount_niche: {min_ncount_niche}. They will be filtered out.")
    niche_assay <- subset(niche_assay, subset = nCount_niche >= 3)
    print(niche_assay)
    write_rds(niche_assay, file.path(dir_res, "niche_sr_object.rds"))
}

if (T) {    
    for (z in c('nCount_niche', 'nFeature_niche', 'nCount_Xenium', 'nFeature_Xenium')) {
        print(z)
        capture.output(quantile(niche_assay@meta.data[, z]), 
        file = file.path(dir_res, glue("quick_stats.{z}.quantiles.txt")))
    }
    p <- VlnPlot(niche_assay, features = c('nCount_niche', 'nFeature_niche', 'nCount_Xenium', 'nFeature_Xenium'), pt.size=0) +
        rremove('x.text') + rremove('x.ticks') + rremove('legend') + rremove('x.title')
    ggsave(file.path(dir_res, "niche_qc_metrics.vlnplot.pdf"), plot=p, width = 6, height = 4, useDingbats = FALSE)
}


#------------------ ~~~ Project to subset ~~~ --------------------
cli_h1("Project to subset")
## The purpose is to reduce the computation burden
niche_assay <- FindVariableFeatures(niche_assay, assay = "niche")
str(VariableFeatures(niche_assay))
df_hvfinfo <- HVFInfo(niche_assay)
table(df_hvfinfo$variance == 0)
feature_variance_zero <- rownames(df_hvfinfo)[which(df_hvfinfo$variance == 0)]
VariableFeatures(niche_assay) <- setdiff(VariableFeatures(niche_assay), feature_variance_zero)
# ncol(niche_assay)
f_niche_sketch_obj <- file.path(dir_res, "niche_sr_object.sketch.rds")
f_niche_ready_obj <- file.path(dir_res, "niche_sr_object.ready.rds")
# if (!file.exists(f_niche_sketch_obj)) {
if (TRUE) {
    n_sketch_cells <- 200e3
    if (ncol(niche_assay) < n_sketch_cells) {
        n_sketch_cells <- round(ncol(niche_assay)/2)
    }
    cli_alert_info("Sketching to {n_sketch_cells} cells")
    tic("SketchData")
    niche_assay <- SketchData(niche_assay, assay = "niche", n = n_sketch_cells)
    toc()
    write_rds(niche_assay, f_niche_sketch_obj)
} else {
    niche_assay <- read_rds(f_niche_sketch_obj)
}
print(niche_assay)
#------------------ ~~~ Preprocessing ~~~ --------------------
cli_h1("Preprocessing")
# DefaultAssay(niche_assay) <- "niche"
DefaultAssay(niche_assay) <- "sketch"
niche_assay <- FindVariableFeatures(niche_assay)
niche_assay <- ScaleData(niche_assay, vars.to.regress = c("nCount_niche"))
cli_alert_info(sprintf("Assay: %s", DefaultAssay(niche_assay)))

if (!'pca' %in% names(niche_assay@reductions)) {
    tic("RunPCA")
    niche_assay <- RunPCA(niche_assay, npcs = 30, verbose = TRUE)
    toc()
    write_rds(niche_assay, f_niche_sketch_obj)
    print(niche_assay)
}
umap_return_model <- FALSE
if (! 'umap' %in% names(niche_assay@reductions) ) {
    if (DefaultAssay(niche_assay) == 'sketch') {
        umap_return_model <- TRUE
    }
    tic("RunUMAP")
    niche_assay <- RunUMAP(niche_assay, dims = 1:30,  return.model = umap_return_model)
    toc()
    write_rds(niche_assay, f_niche_sketch_obj)
    print(niche_assay)
}

#------ Quick check if clusters are driven by samples ------
if (T) {
    pal_samples <- init_pal_d(niche_assay@meta.data$sample, pal = "parade")
}
for (z in c('sample', 'celltype')) {
    pal_use <- switch(z,
        sample = pal_samples,
        celltype = pal_celltypes)
    p <- UMAPPlot(niche_assay, group.by = z, cols = pal_use, label=TRUE, pt.size=3, raster=TRUE, raster.dpi = c(1024, 1024)) +
        ggtitle(z) +
        NoLegend()
    ggsave(file.path(dir_res, glue("umap.{z}.pdf")), plot=p, width = 7, height = 7, useDingbats = FALSE)
}
for (z in c('nCount_niche', 'nFeature_niche', 'nCount_Xenium', 'nFeature_Xenium')) {
    p <- FeaturePlot(niche_assay, features = z, pt.size=3, raster=TRUE, raster.dpi = c(1024, 1024), order = TRUE) +
        ggtitle(z)
    ggsave(file.path(dir_res, glue("umap.{z}.pdf")), plot=p, width = 7, height = 7, useDingbats = FALSE)
}
#------------------ ~~~ Clustering ~~~ --------------------
cli_h1("Clustering")
if (is.null(names(niche_assay@graphs))) {
    tic("FindNeighbors")
    niche_assay <- FindNeighbors(niche_assay, reduction = "pca", dims = 1:30)
    toc()
    write_rds(niche_assay, f_niche_sketch_obj)
    print(niche_assay)
}

cat("[done]")
timestamp()
