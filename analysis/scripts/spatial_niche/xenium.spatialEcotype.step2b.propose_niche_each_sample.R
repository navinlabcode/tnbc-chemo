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
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.nmf.viz.R")
    library(openxlsx)
    options(ggrastr.default.dpi = 600)
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    f_assay <- cmdargs[1]
    f_xenium <- cmdargs[2] # for getting sample name
} else {
    sample <- "ART235"
    f_assay <- glue("/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/ALL/scNicheRadius/R30/split_sample/{sample}/niche_proposal_onestep/niche_sr_object.sketch.rds")
    f_xenium <- glue("/volumes/USR1/yyan/project/tnbc_xenium/data/{sample}/cleaned1/ready.seurat.rds")
}

#------ Libs ------
if (T) {
    npcs <- 30
}

f_dict_snn_to_niche <- "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/ALL/scNicheRadius/R30/split_sample/dataframe.each_sample_niche_k.exported.csv"
    sample_array <- c(
        "ART10", "ART23", "ART312", "ART3122", "ART311", "ART305", "ART304",
        "ART18", "ART31", "ART40", "ART43", "ART65",
        "ART94", "ART133",
        "ART117", "ART217", "ART258", "ART272", "ART232", "ART282", "ART289", "ART122", "ART250", "ART92",
        "ART153", "ART194", "ART238", "ART257", "ART104", "ART279", "ART30", "ART71", "ART247", "ART271",
        "ART155", "ART170", "ART219", "ART236", "ART277", "ART202", "ART223", "ART235", "ART266", "ART276"
    )

if (! file.exists(f_dict_snn_to_niche)) {
    dict_sample_niche_k <- data.frame(
        sample = sample_array,
        niche_k = 8
    ) %>% arrange(sample) ## initial a worksheet to work later
    openxlsx::write.xlsx(
        dict_sample_niche_k,
        file.path(
            "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44/spatial_ecotype_winner/ALL/scNicheRadius/R30/split_sample",
            "dataframe.each_sample_niche_k.xlsx"
        )
    )
} else {
    dict_sample_niche_k <- read_csv(f_dict_snn_to_niche, show_col_types = FALSE) %>%
        arrange(sample)
    print(dict_sample_niche_k)
}

#------------------ ~~~ Inputs ~~~ --------------------
cli_h1("Inputs")

guess_sample <- str_extract(f_assay, "ART[0-9]+")

if (! guess_sample %in% sample_array) {
    stop(glue("Sample {guess_sample} not in the allowed sample_array!"))
} else {
    sample <- guess_sample
    message(glue("Sample: {sample}"))
}

obj <- readRDS(f_assay)
print(obj)
print(DefaultAssay(obj))

xmo <- readRDS(f_xenium)
print(xmo)
print(DefaultAssay(xmo))

dir_assay <- dirname(f_assay); dir_res <- file.path(dir_assay, "snn_clusters")
fs::dir_create(dir_res)

#------ Spatial preparation ------
if (T) {
    img_pt_size <- .2
    insitu_draw_dark_bg <- FALSE
    xmo_coord <- GetTissueCoordinates(xmo)
    tissue_xmax <- max(xmo_coord$x)
    tissue_xmin <- min(xmo_coord$x)
    tissue_ymax <- max(xmo_coord$y)
    tissue_ymin <- min(xmo_coord$y)
    tissue_width <- tissue_xmax - tissue_xmin
    tissue_height <- tissue_ymax - tissue_ymin
    tissue_xyratio <- tissue_width / tissue_height
    tissue_img_width <- calc_img_spatial_pdf_size(tissue_xyratio)["width"]
    tissue_img_height <- calc_img_spatial_pdf_size(tissue_xyratio)["height"]
    # tissue_img_width <- to_img_xenium_pdf_inch(tissue_width)
    # tissue_img_height <- to_img_xenium_pdf_inch(tissue_height)
    ncells_tissue <- sum(xmo_coord$x >= tissue_xmin & xmo_coord$x <= tissue_xmax &
        xmo_coord$y >= tissue_ymin & xmo_coord$y <= tissue_ymax)

    magnitude_ROI_size <- log10(max(tissue_width, tissue_height))

    img_pt_cex <- case_when(
        ncells_tissue > 0 ~ 4,
        ncells_tissue > 20 ~ 2,
        ncells_tissue > 30 ~ 1.5,
        ncells_tissue > 5000 ~ 1,
        ncells_tissue > 10000 ~ 0.5,
        ncells_tissue > 20000 ~ 0.35,
        ncells_tissue > 50000 ~ 0.2,
        .default = 0
    )
    # img_pt_cex <- min(tissue_width, tissue_height) / 1000 / PROBE_PT_CEX_DENOMINATOR
    img_pt_cex <- 1 + 2^(-1 * magnitude_ROI_size)

    probe_pt_size <- img_pt_size / 10
    probe_pt_cex <- img_pt_cex # min(tissue_width, tissue_height) / 1000 / PROBE_PT_CEX_DENOMINATOR
}

#------------------ ~~~ Clustering on sketch data ~~~ --------------------
cli_h1("Clustering on sketch data")
res_tried <- c(0.2, 0.5, 1, 2)
for (res in res_tried) {
    message(glue::glue("Resolution: {res}"))
    obj <- FindClusters(obj, resolution = res)
    res_str <- paste0("sketch_snn_res.", res)
    obj@meta.data[, res_str] <- obj@meta.data[, "seurat_clusters"] ## pretify the factor order
}
res_sketch_string <- paste0("sketch_snn_res.", res_tried)
all(res_sketch_string %in% colnames(obj@meta.data))
colnames(obj@meta.data)
table(obj$seurat_clusters, useNA = "ifany")

#------------------ ~~~ Project to the full dataset ~~~ --------------------
cli_h1("Project to the full dataset")
for (res_str in res_sketch_string) {
    message(glue::glue("Projecting Resolution: {res_str}"))
    print(table(Current = obj@meta.data[, res_str], useNA = "ifany"))
    if (!any(is.na(obj@meta.data[, res_str]))) {
        next()
    }
    obj <- ProjectData(
        object = obj,
        assay = "niche",
        full.reduction = "pca.full",
        sketched.assay = "sketch",
        sketched.reduction = "pca",
        umap.model = "umap",
        dims = 1:npcs,
        refdata = list(cluster_full = res_str)
    )
    obj@meta.data[, res_str] <- obj@meta.data[, "cluster_full"]
    obj@meta.data[, "cluster_full"] <- NULL
    print(table(After = obj@meta.data[, res_str], useNA = "ifany"))
}

for (res_str in res_sketch_string) {
    ## prefify the factor order
    obj@meta.data[, res_str] <- factor(
        as.character(obj@meta.data[, res_str]),
        levels = gtools::mixedsort(unique(as.character(obj@meta.data[, res_str])))
    )
}
write_rds(obj@meta.data[, res_sketch_string], file.path(dir_res, "dataframe.cellmeta_niche_sketch_snn.rds"))

#------ OPTIONS: other clustering methods ------
# Rphenograph
# Simple k-means or Hierarchical clustering
# Density-Based Spatial Clustering of Applications with Noise (DBSCAN)

write_rds(obj@meta.data, file.path(dir_res, "dataframe.cellmeta.allcolumns.rds"))


#------ Legend ------
pal_list <- list()
for (res_str in res_sketch_string) {
    pal_res_str <- init_pal_d(obj@meta.data[, res_str])
    pal_res_str["NA"] <- "ghostwhite"
    pal_to_ggplot(pal_res_str, pal_name = res_str) %>% ggsave(
        filename = file.path(dir_res, glue("legend_{res_str}.pdf")),
        plot = ., width = 2, height = length(pal_res_str) * 0.25 + 0.5, useDingbats = FALSE
    )
    pal_list[[res_str]] <- pal_res_str
}

#------ Dimplot of niche clusters ------
cli_h1("DimPlot of sketch clusters")
for (res_str in res_sketch_string) {
    cli_h2(res_str)
    p <- DimPlot(
        obj,
        group.by = res_str, cols = pal_list[[res_str]],
        label = TRUE, label.size = 1, repel = TRUE, raster = TRUE,
        pt.size = 2, raster.dpi = c(1024, 1024)
    ) +
        theme_void() +
        theme(aspect.ratio = 1) +
        labs(title = glue("res={res_str}")) +
        rremove("legend")
    ggsave(
        filename = file.path(dir_res, glue("umap_sketch_snn_clusters_{res_str}.pdf")),
        plot = p,
        width = 3.5, height = 3, useDingbats = FALSE
    )
}

#------ SpatialDimplt of niche clusters ------
cli_h1("SpatialDimPlot of sketch clusters")

for (res_str in res_sketch_string) {
    cli_h2(res_str)
    pal_res_str <- pal_list[[res_str]]
    xmo <- AddMetaData(xmo, metadata = obj@meta.data[colnames(xmo), res_str, drop = TRUE], col.name = res_str)

    # p <- ImageDimPlot.raster(
    #     xmo,
    #     group.by = res_str, cols = pal_res_str,
    #     size = img_pt_size, cex = img_pt_cex) +
    #     labs(title = glue("Sketch SNN clusters (res={res_str})"))

    p <- ImageDimPlotSegmentation.raster(
        xmo,
        group.by = res_str, cols = pal_res_str,
        cell_segmentation_width = 0.01
    ) +
        labs(title = glue("res={res_str}")) +
        rremove("legend")
    p <- q_add_anno_xenium_scale_bar(
        p,
        roi_xmin = tissue_xmin, roi_xmax = tissue_xmax, roi_ymin = tissue_ymin, roi_ymax = tissue_ymax,
        color = ifelse(insitu_draw_dark_bg, "white", "black"), linewidth = 3
    )

    ggsave(
        filename = file.path(dir_res, glue("spatial_sketch_snn_clusters_{res_str}.pdf")),
        plot = p,
        width = tissue_img_width, height = tissue_img_height, useDingbats = FALSE
    )
}

#------------------ ~~~ Annotate niche clusters ~~~ --------------------
cli_h1("Annotate niche clusters")
BIO_CELL_IDENT_STR <- "cell_state_paper"
freq_BIO_CELL_IDENT <- c(prop.table(table(obj@meta.data[, BIO_CELL_IDENT_STR])))
names(freq_BIO_CELL_IDENT) <- str_remove_all(names(freq_BIO_CELL_IDENT), "module_score_")
res_str <- "sketch_snn_res.2"
table(obj@meta.data[[res_str]], useNA = "ifany")
prop.table(table(obj@meta.data[[res_str]], useNA = "ifany"))

obj_avg <- AverageExpression(
    obj,
    assays = "niche",
    features = rownames(obj),
    group.by = res_str,
    slot = "counts",
    return.seurat = TRUE
)
obj_avg <- NormalizeData(obj_avg, normalization.method = "RC", scale.factor = 100)
obj_avg <- ScaleData(obj_avg, features = rownames(obj_avg), scale.max = 4)
# obj_avg <- ScaleData(obj_avg, features = rownames(obj_avg))
write_rds(obj_avg, file.path(dir_res, glue("niche_avg_object_{res_str}.rds")))


#------ Merge niche clusters into niche ------
### Generate a dictionary to map res_str to niche
### Methods: ConsensusClusterPlus

# obj_avg_mat <- obj_avg_mat_freq_logenrich
obj_avg_mat <- as.matrix(GetAssayData(obj_avg, assay = "niche", slot = "scale.data"))
write_rds(obj_avg_mat, file.path(dir_res, glue("niche_avg_scale_matrix_{res_str}.rds")))
# obj_avg_mat <- readRDS(file.path(dir_res, glue("niche_avg_scale_matrix_{res_str}.rds")))
# obj_avg_mat <- as.matrix(GetAssayData(obj_avg, assay = "niche", slot = "data"))

sum(obj_avg_mat[, 2], na.rm = TRUE)
range(obj_avg_mat)
n_cluster <- ncol(obj_avg_mat)
n_feature <- nrow(obj_avg_mat)

library(ConsensusClusterPlus)
# Run ConsensusClusterPlus to determine the optimal number of clusters
dir_CSSCluster <- file.path(dir_res, glue("ConsensusClusterPlus.{res_str}"))
fs::dir_create(dir_CSSCluster)
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.ConsensusClusterPlus.R")
maxK <- 15
css_cluster_results <- ConsensusClusterPlus(
    d = obj_avg_mat,
    maxK = maxK, # Maximum number of clusters to evaluate
    reps = 100, # Number of resampling iterations
    pItem = 0.8, # Proportion of items to sample
    pFeature = 1, # Proportion of features to sample
    clusterAlg = "hc", # Clustering algorithm (hierarchical clustering)
    distance = "euclidean", # Distance metric
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    corUse = "pairwise.complete.obs",
    seed = 42, # Set seed for reproducibility
    plot = "pdf", # Output plots as PDF
    title = dir_CSSCluster # Output directory for result
)
write_rds(css_cluster_results, file.path(dir_CSSCluster, "ConsensusClusterPlus_results.rds"))
# str(css_cluster_results)
# length(css_cluster_results)
if (TRUE) {
    CDF_values <- readRDS(file.path(dir_CSSCluster, "CDF_values.rds"))
    CDF_cumsum <- cumsum(CDF_values)
    ## find the elbow point
    library(pathviewr)
    elbow_point_xx <- ElbowPoint_gives_optK(k = c(1:length(CDF_cumsum)) + 1, v = CDF_cumsum)
    print(elbow_point_xx)
    pdf(file.path(dir_CSSCluster, "CDF_cumsum_elbow.pdf"), width = 4, height = 4)
    plot((1:length(CDF_cumsum)) + 1, CDF_cumsum,
        type = "b",
        xlab = "Number of clusters (K)", ylab = "Cumulative change in area under CDF curve",
        main = glue("Delta Area best={elbow_point_xx}")
    )
    abline(v = elbow_point_xx, col = "red", lty = 2)
    dev.off()
}

library(factoextra)
pdf(file.path(dir_CSSCluster, "NbClust_GapStat.pdf"), width = 4, height = 4)
p <- fviz_nbclust(x = t(obj_avg_mat), FUNcluster = hcut, method = "gap_stat", k.max = maxK)
# p <- fviz_nbclust(x = t(obj_avg_mat), FUNcluster = hcut, method = "silhouette", k.max = maxK)
print(p)
dev.off()

# geom_vline(xintercept = elbow_point_xx, linetype = 2, color = "red") +
# labs(title = glue("Elbow method best={elbow_point_xx}"))

#------ Heatmap at its own ------
ht <- ComplexHeatmap::Heatmap(
    obj_avg_mat,
    name = "zscore",
    column_split = 8,
    cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
    cluster_columns = TRUE, show_column_dend = TRUE, column_names_side = "top", column_dend_height = unit(5, "mm"),
    column_title = res_str,
    column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
    # col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    col = heatmap_color_fun_zscore_2,
    # col = heatmap_color_fun_cont(from=quantile(obj_avg_mat, 0.01), to=quantile(obj_avg_mat, 0.99), palette='Blue-Red'),
    # col = viridis::viridis(256),
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
    width = unit(5 * n_cluster / n_feature, "inches"),
    height = unit(5, "inches"), use_raster = F
)
pdf(file.path(dir_res, glue("heatmap_niche_avg_{res_str}.pdf")),
    width = 5 * n_cluster / n_feature + 1, height = 5 + 2,
    useDingbats = FALSE
)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

#------ Heatmap with column order by CSSCluster ------
str(css_cluster_results[[2]])

df_css_cluster <- lapply(2:maxK, function(k) {
    o <- data.frame(
        sample = names(css_cluster_results[[k]]$consensusClass),
        res = css_cluster_results[[k]]$consensusClass
    )
    colnames(o)[2] <- paste0("cut", k)
    o
}) %>% Reduce(left_join, .)
# view(df_css_cluster)
column_to_rownames(df_css_cluster, "sample") -> df_css_cluster
write_rds(df_css_cluster, file.path(dir_res, glue("dataframe.cellmeta_css_clusters_{res_str}.rds")))

pal_css_cluster <- lapply(2:maxK, function(k) {
    cl <- css_cluster_results[[k]]$clrs[[3]]
    if (length(cl) < k) {
        cl <- rainbow(n=k)
    }
    o <- structure(cl, names = 1:length(cl))
    o
})
names(pal_css_cluster) <- paste0("cut", 2:maxK)
write_rds(pal_css_cluster, file.path(dir_res, glue("palette_css_clusters_{res_str}.rds")))
# pal_css_cluster <- readRDS(file.path(dir_res, glue("palette_css_clusters_{res_str}.rds")))

ht_col_order <- css_cluster_results[[maxK]]$consensusTree$order

optimal_cut <- dict_sample_niche_k[dict_sample_niche_k$sample == sample, "niche_k", drop=T]
cli_h2(glue("Optimal cut: {optimal_cut}"))
ComplexHeatmap::HeatmapAnnotation(
    df = df_css_cluster[ht_col_order, ],
    col = pal_css_cluster,
    which = "column",
    show_legend = F, annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size = unit(0.08, "inches"), gap = 0
    # annotation_height = unit(1, "inches"),
    # height = unit(9, "inches")
) -> ha_top
ht_css <- ComplexHeatmap::Heatmap(
    obj_avg_mat[, ht_col_order],
    name = "zscore",
    top_annotation = ha_top,
    cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
    cluster_columns = TRUE, show_column_dend = FALSE, column_names_side = "top", column_dend_height = unit(5, "mm"),
    column_split = optimal_cut, column_gap = unit(0, "mm"), border = TRUE,
    column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
    column_title = res_str,
    col = heatmap_color_fun_zscore_2,
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
    width = unit(5 * n_cluster / n_feature, "inches"),
    height = unit(5, "inches"), use_raster = F
)
pdf(file.path(dir_res, glue("heatmap_css_niche_avg_{res_str}.pdf")),
    width = 5 * n_cluster / n_feature + 1, height = 5 + 1 + 0.08 * maxK,
    useDingbats = FALSE
)
draw(ht_css, heatmap_legend_side = "bottom")
dev.off()

#------------------ ~~~ Determine the best cut ~~~ --------------------
cli_h1("Determine the best k")
optimal_cut <- dict_sample_niche_k[dict_sample_niche_k$sample == sample, "niche_k", drop=T]
dict_snn_to_css_optimal <- structure(
    df_css_cluster[[paste0("cut", optimal_cut)]], 
    names = rownames(df_css_cluster))
dict_snn_to_css_optimal <- df_css_cluster[, paste0("cut", optimal_cut), drop=F] %>% rownames_to_column("snn")
colnames(dict_snn_to_css_optimal) <- c("snn", "niche_cluster")

obj@meta.data[, 'snn'] <- paste0('g', as.character(obj@meta.data[, res_str]))
cellmeta_niche_cluster <- left_join(
    obj@meta.data[, 'snn', drop=F], 
    dict_snn_to_css_optimal
)
obj@meta.data[, 'niche_cluster'] <- cellmeta_niche_cluster$niche_cluster

write_rds(obj@meta.data, file.path(dir_res, "dataframe.cellmeta.allcolumns.with_niche_cluster.rds"))

#------ Heatmap showing the niche cluster (signature) ------
pal_niche <- pal_css_cluster[[paste0("cut", optimal_cut)]]
write_rds(pal_niche, file.path(dir_res, glue("palette_niche_clusters_{optimal_cut}.rds")))
ggsave(
    pal_to_ggplot(pal_niche, pal_name = glue("niche_cluster_k{optimal_cut}")),
    filename = file.path(dir_res, glue("legend_niche_cluster_k{optimal_cut}.pdf")),
    width = 2, height = length(pal_niche) * 0.25 + 0.5, useDingbats = FALSE
)

niche_str <- "niche_cluster"
obj_niche <- AverageExpression(
    obj,
    assays = "niche",
    features = rownames(obj),
    group.by = niche_str,
    slot = "counts",
    return.seurat = TRUE
) ## To clarify, obj_niche's counts are not necessary integers

if (F) {
    ## Check if the AverageExpression is done as my expectation [Yes]
    yy <- as.matrix(GetAssayData(obj_niche, assay = "niche", slot = "counts"))
    yy[1:3, 1:3]
    dim(yy)
    xx <- ruok::mat_tapply(
        mat=t(GetAssayData(obj, assay = "niche", slot = "counts")), 
        INDEX = paste0('g', obj@meta.data[[niche_str]]), 
        FUN = mean, rm.na=TRUE)
    xx <- t(xx)
    xx[1:3, 1:3]
    range(xx)
    range(yy)
    all.equal(yy, xx)
}

obj_niche <- NormalizeData(obj_niche, normalization.method = "RC", scale.factor = 100)
obj_niche <- ScaleData(obj_niche, features = rownames(obj_niche), scale.max = 4)
write_rds(obj_niche, file.path(dir_res, glue("niche_avg_object_{niche_str}.rds")))
obj_niche_mat <- as.matrix(GetAssayData(obj_niche, assay = "niche", slot = "scale.data"))
rownames(obj_niche_mat) <- str_remove_all(rownames(obj_niche_mat), "module_score_")
write_rds(obj_niche_mat, file.path(dir_res, glue("niche_avg_scale_matrix_{niche_str}.rds")))
# obj_niche_mat <- readRDS(file.path(dir_res, glue("niche_avg_scale_matrix_{niche_str}.rds")))
dim(obj_niche_mat)
colnames(obj_niche_mat) <- str_replace_all(colnames(obj_niche_mat), replacement = "n", "g")

ht_col_anno <- ComplexHeatmap::HeatmapAnnotation(
    niche = colnames(obj_niche_mat),
    col = list(niche = structure(as.character(pal_niche), names = paste0('n', names(pal_niche)))),
    which = "column", show_annotation_name = F,
    show_legend = F, annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size = unit(2, "mm"), gap = 0
    # annotation_height = unit(1, "inches"),
    # height = unit(9, "inches")
)
ht <- ComplexHeatmap::Heatmap(
    obj_niche_mat,
    name = "zscore",
    top_annotation = ht_col_anno,
    cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
    cluster_columns = TRUE, show_column_dend = TRUE, column_names_side = "top", column_dend_height = unit(3, "mm"),
    column_title = niche_str,
    column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
    col = heatmap_color_fun_zscore_2,
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
    width = unit(5 * ncol(obj_niche_mat) / nrow(obj_niche_mat), "inches"),
    height = unit(5, "inches"), use_raster = F
)
pdf(file.path(dir_res, glue("heatmap_niche_avg_{niche_str}.pdf")),
    width = 5 * ncol(obj_niche_mat) / nrow(obj_niche_mat) + 1, height = 5 + 2,
    useDingbats = FALSE
)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

#------ Dimplot ------
p <- DimPlot(
    obj,
    group.by = niche_str, cols = pal_niche,
    label = TRUE, label.size = 3, repel = TRUE, raster = TRUE,
    pt.size = 2, raster.dpi = c(1024, 1024)
) +
    theme_void() +
    theme(aspect.ratio = 1) +
    labs(title = glue("Niche clusters (k={optimal_cut})")) +
    rremove("legend")
ggsave(
    filename = file.path(dir_res, glue("umap_niche_clusters_{niche_str}.pdf")),
    plot = p,
    width = 3.5, height = 3, useDingbats = FALSE
)
#------ SpatialDimPlot ------
pal_niche["NA"] <- "ghostwhite"
xmo <- AddMetaData(xmo, metadata = obj@meta.data[colnames(xmo), niche_str, drop = TRUE], col.name = niche_str)
p <- ImageDimPlotSegmentation.raster(
    xmo,
    group.by = niche_str, cols = pal_niche,
    cell_segmentation_width = 0.01
) +
    labs(title = glue("Niche clusters (k={optimal_cut})")) +
    rremove("legend")
p <- q_add_anno_xenium_scale_bar(
    p,
    roi_xmin = tissue_xmin, roi_xmax = tissue_xmax, roi_ymin = tissue_ymin, roi_ymax = tissue_ymax,
    color = ifelse(insitu_draw_dark_bg, "white", "black"), linewidth = 3
)
ggsave(
    filename = file.path(dir_res, glue("spatial_niche_clusters_{niche_str}.pdf")),
    plot = p,
    width = tissue_img_width, height = tissue_img_height, useDingbats = FALSE
)



cat("All done!\n")
