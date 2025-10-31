suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: QC and prepare the raw Seurat object for Visium HD data
#       and save it as a RDS file.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025-05-13
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
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.visiumHD.R")
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    # f_in <- cmdargs[1]
    sample_name <- cmdargs[1]
    vs_binsize <- as.numeric(cmdargs[2])
} else {
    f_in <- "/volumes/USR1/yyan/project/tnbc_visium_hd/data0/ART122/pass.seurat.rds"
    sample_name <- "ART122"
    vs_binsize <- 32
}
f_in <- file.path(
    "/volumes/USR1/yyan/project/tnbc_visium_hd/data0",
    sample_name, sprintf("binsize_%s", vs_binsize),
    "pass.seurat.rds"
)
dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_visium_hd/data", sample_name,
    sprintf("binsize_%s", vs_binsize)
)
fs::dir_create(dir_res)

param.min.nCount <- 10
param.min.nFeature <- 10
if (vs_binsize == 32) {
    param.min.nCount <- 10
    param.min.nFeature <- 30
}
if (vs_binsize == 8) {
    param.min.nCount <- 10
    param.min.nFeature <- 10
}
#------------------ ~~~ Read in ~~~ --------------------
cli_h1("Read in")

sr3 <- read_rds(f_in)
print(sr3)
DefaultAssay(sr3)
print(range(sr3$nFeature_Spatial))
print(range(sr3$nCount_Spatial))

fail_QC <- !(sr3$nFeature_Spatial >= param.min.nFeature & sr3$nCount_Spatial >= param.min.nCount)
fail_QC <- ifelse(fail_QC, "Fail", "Pass")
sr3$fail_QC <- fail_QC
print(table(sr3$fail_QC))

sp_xy_ratio <- get_visium_xy_ratio(sr3)
pdf_width <- 5
pdf_height <- 5
if (sp_xy_ratio > 1) {
    pdf_height <- pdf_height * sp_xy_ratio
} else {
    pdf_width <- pdf_width * sp_xy_ratio
}
print(c(pdf_width, pdf_height))


p <- SpatialDimPlot(
    sr3,
    group.by = "fail_QC", pt.size.factor = 2,
    cols = c("Fail" = "red", "Pass" = "blue")
) + theme(aspect.ratio = sp_xy_ratio) +
    labs(subtitle = sprintf(
        "%s:  %s/%s cells are kept",
        sample_name, comma(sum(sr3$fail_QC == "Pass")), comma(length(sr3$fail_QC))
    )) +
    theme(legend.position = "none") +
    ggtitle("QC of cells")
ggsave(
    filename = file.path(dir_res, "spatial.QC_cells.pdf"),
    plot = p,
    width = pdf_width,
    height = pdf_height, useDingbats = F
)

#------------------ ~~~ Subset ~~~ --------------------
cli_h1("Subset")
sr3 <- subset(sr3, subset = fail_QC == "Pass")
print(sr3)

#------------------ ~~~ Prepare ~~~ --------------------
cli_h1("Prepare")
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")

f_out <- file.path(dir_res, "ready.seurat.rds")
if (!file.exists(f_out)) {
    cli_alert_info("File {f_out} does not exist, start the preparation.")

    sr3 <- FindVariableFeatures(sr3, nfeatures = 3000)
    sr3 <- NormalizeData(sr3)
    sr3 <- ScaleData(sr3)
    sr3 <- RunPCA(sr3, reduction.name = "pca", npcs = 50)
    sr3 <- FindNeighbors(sr3, reduction = "pca", dims = 1:50)

    sr3 <- FindClusters(sr3, resolution = 0.2)
    sr3 <- RunUMAP(sr3, reduction = "pca", reduction.name = "umap", dims = 1:50)
    write_seurat(sr3, dir_res, "ready")
} else {
    cli_alert_info("File {f_out} exists, skip the preparation.")
    sr3 <- read_rds(f_out)
}
# Spatial_snn_res.0.6
print(colnames(sr3@meta.data))

#------ export the UMI matrix ------
f_umi <- file.path(dir_res, "umi_count.matrix.rds")
if (!file.exists(f_umi)) {
    cli_alert_info("Save the UMI count matrix")
    sr3 %>%
        GetAssayData(slot = "counts") %>%
        readr::write_rds(f_umi)
}
#------ Clusters ------
cli_h2("Plotting clusters")
z <- "Spatial_snn_res.0.2"
p <- SpatialDimPlot(
    sr3,
    group.by = z, pt.size.factor = 5,
    label = F, alpha = 1,
) + theme(aspect.ratio = sp_xy_ratio)
ggsave(
    filename = file.path(dir_res, sprintf("spatial.%s.pdf", z)),
    plot = p,
    width = pdf_width + 3,
    height = pdf_height, useDingbats = F
)
p <- DimPlot(
    sr3,
    group.by = z,
    label = T, alpha = 1, raster = T, pt.size = 3
) + theme(aspect.ratio = 1)
ggsave(
    filename = file.path(dir_res, sprintf("dimplot.%s.pdf", z)),
    plot = p,
    width = 4,
    height = 4, useDingbats = F
)
#------ Genes ------
cli_h2("Plotting genes")
g_sets <- c(
    "EPCAM", "PTPRC",
    "CD3D", "BANK1", "NKG7",
    "CD68", "FCGR3A", "CD163",
    "LUM", "VWF"
)
fs::dir_create(file.path(dir_res, "specific_genes"))
library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = brewer.pal(n = 9, name = "YlOrRd"))

for (g in g_sets) {
    p <- SpatialFeaturePlot(
        sr3,
        features = g, pt.size.factor = 5,
        min.cutoff = "q1", max.cutoff = "q99",
        alpha = 1
    ) + theme(aspect.ratio = sp_xy_ratio)
    # scale_fill_gradientn(colors = SpatialColors(100))
    ggsave(
        filename = file.path(dir_res, "specific_genes", sprintf("spatial.feature.%s.pdf", g)),
        plot = p,
        width = pdf_width + 3,
        height = pdf_height, useDingbats = F
    )
    p <- FeaturePlot(
        sr3,
        features = g,
        min.cutoff = "q1", max.cutoff = "q99",
        alpha = 1, raster = T, pt.size = 3
    ) + theme(aspect.ratio = 1)
    # scale_color_gradientn(colors = SpatialColors(100))

    ggsave(
        filename = file.path(dir_res, "specific_genes", sprintf("dimplot.feature.%s.pdf", g)),
        plot = p,
        width = 4,
        height = 4, useDingbats = F
    )
}
#------ nFeatures and nCounts ------
cli_h2("Plotting nFeature and nCount")
sr3$log2nCount_Spatial <- log2(sr3$nCount_Spatial + 1)
for (y in c("nFeature_Spatial", "nCount_Spatial", "log2nCount_Spatial")) {
    p <- SpatialFeaturePlot(
        sr3,
        features = y, pt.size.factor = 5,
        min.cutoff = "q1", max.cutoff = "q99",
        alpha = 1
    ) + theme(aspect.ratio = sp_xy_ratio)
    ggsave(
        filename = file.path(dir_res, sprintf("spatial.feature.%s.pdf", y)),
        plot = p,
        width = pdf_width + 3,
        height = pdf_height, useDingbats = F
    )
    p <- FeaturePlot(
        sr3,
        features = y,
        min.cutoff = "q1", max.cutoff = "q99",
        label = T, alpha = 1, raster = T, pt.size = 3
    ) + theme(aspect.ratio = 1)
    ggsave(
        filename = file.path(dir_res, sprintf("dimplot.feature.%s.pdf", y)),
        plot = p,
        width = 4,
        height = 4, useDingbats = F
    )
}
#------ Just the image ------
p <- SpatialDimPlot(
    sr3,
    group.by = "orig.ident", pt.size.factor = 1,
    label = F, alpha = 0
) + theme(aspect.ratio = sp_xy_ratio) +
    labs(subtitle = sprintf("%s", sample_name)) +
    theme(legend.position = "none")
ggsave(
    filename = file.path(dir_res, "spatial.image.pdf"),
    plot = p,
    width = pdf_width + 3,
    height = pdf_height, useDingbats = F
)

cat("[done]")
timestamp()
