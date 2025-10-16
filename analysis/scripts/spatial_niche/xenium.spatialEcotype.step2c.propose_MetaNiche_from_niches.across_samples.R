suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal:
# Concatenate the SNN clusters of all samples to find the similar clusters across samples.
# These similar clusters are proposed as Meta Niches.
#
# Input: the niche score matrix of SNN clusters in each sample (output from xenium.spatialEcotype.step2b.propose_niche_each_sample.R)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025-08-29
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
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.nmf.viz.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
    options(ggrastr.default.dpi = 600)
    library(Seurat)
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    dir_proj <- cmdargs[1]
} else {
    dir_proj <- file.path(
        "/volumes/USR1/yyan/project/tnbc_xenium/data_merged_N44",
        "spatial_ecotype_winner", ## change here for different methods of assigning cell labels
        "ALL", ## change here for different contexts; normally do not change
        "scNicheRadius",
        "R30" ## change here for different radius; normally do not change
    )
}
dir_res <- file.path(dir_proj, "MetaNiche_across_samples")
fs::dir_create(dir_res)

sample_array <- c(
    "ART10", "ART23", "ART312", "ART3122", "ART311", "ART305", "ART304",
    "ART18", "ART31", "ART40", "ART43", "ART65",
    "ART94", "ART133",
    "ART117", "ART217", "ART258", "ART272", "ART232", "ART282", "ART289", "ART122", "ART250", "ART92",
    "ART153", "ART194", "ART238", "ART257", "ART104", "ART279", "ART30", "ART71", "ART247", "ART271",
    "ART155", "ART170", "ART219", "ART236", "ART277", "ART202", "ART223", "ART235", "ART266", "ART276"
)
str(unique(sample_array))

#------------------ ~~~ Inputs ~~~ --------------------
cli_h1("Inputs")
res_str <- "niche_cluster"

f_mat_score <- file.path(dir_res, glue("niche_avg_scale_matrix_{res_str}_all_samples.rds"))
f_mat_data <- file.path(dir_res, glue("niche_avg_data_matrix_{res_str}_all_samples.rds"))

# if (!file_exists(f_mat_score)) {
if (T) {
    file.path(
        dir_proj,
        "split_sample",
        sample_array, "niche_proposal_onestep/snn_clusters/",
        glue("niche_avg_object_{res_str}.rds")
    ) -> file_score_list
    names(file_score_list) <- sample_array
    if (!all(file_exists(file_score_list))) {
        cli_alert_danger("Some input files do not exist!")
        cli_ol(file_score_list[!file_exists(file_score_list)])
    } else {
        cli_alert_success("All input files exist.")
    }

    mat_score_list <- lapply(1:length(file_score_list), function(i) {
        x <- readRDS(file_score_list[i])
        x <- GetAssayData(x, assay = "niche", layer = "scale.data")
        ## each mat is a in a shape of features (cell states) x clusters
        sample_name <- names(file_score_list)[i]
        colnames(x) <- paste0(sample_name, "_", colnames(x))
        return(x)
    })
    print(mat_score_list[[1]][1:2, 1:2])
    union_features <- gtools::mixedsort(unique(unlist(lapply(mat_score_list, rownames))))
    str(union_features)

    mat_score_list <- lapply(mat_score_list, function(x) {
        idx <- match(union_features, rownames(x))
        x <- x[idx, , drop = FALSE]
        rownames(x) <- union_features
        # x[is.na(x)] <- 0
        return(x)
    })

    mat_score <- do.call(cbind, mat_score_list)
    write_rds(mat_score, f_mat_score)

    assay_list <- lapply(1:length(file_score_list), function(i) {
        x <- readRDS(file_score_list[i])
        s <- names(file_score_list)[i]
        x <- RenameCells(x, new.names = paste0(s, "_", Cells(x)))
        x <- GetAssayData(x, assay = "niche", layer = "data") ## Relative abundance
        x <- as.matrix(x)
        return(x)
    })
    assay_list <- lapply(assay_list, function(x) {
        idx <- match(union_features, rownames(x))
        x <- x[idx, , drop = FALSE]
        rownames(x) <- union_features
        x[is.na(x)] <- 0 ## Relative abundance should be zero
        return(x)
    })
    mat_data <- assay_merged <- do.call(cbind, assay_list)
    write_rds(assay_merged, f_mat_data)
} else {
    mat_score <- readRDS(f_mat_score)
    mat_data <- readRDS(f_mat_data)
}
print(dim(mat_score))
print(mat_score[1:2, 1:2])

#------ Heatmap on its own ------
n_cluster <- ncol(mat_score)
n_feature <- nrow(mat_score)
ht <- ComplexHeatmap::Heatmap(
    mat_score,
    name = "zscore",
    column_split = 8,
    cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
    cluster_columns = TRUE, show_column_dend = TRUE, column_names_side = "top", column_dend_height = unit(5, "mm"),
    column_title = res_str, show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
    col = heatmap_color_fun_zscore_2, na_col = "ghostwhite",
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
    width = unit(7, "inches"),
    height = unit(5, "inches"), raster_by_magick = TRUE
)
pdf(file.path(dir_res, glue("heatmap_niche_score_matrix_all_samples.pdf")),
    width = 7 + 1, height = 5 + 1, useDingbats = FALSE
)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

#------------------ ~~~ Clustering ~~~ --------------------
cli_h1("Clustering")

library(ConsensusClusterPlus)
# Run ConsensusClusterPlus to determine the optimal number of clusters
dir_CSSCluster <- file.path(dir_res, glue("ConsensusClusterPlus.{res_str}"))
fs::dir_create(dir_CSSCluster)
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.ConsensusClusterPlus.R")
maxK <- 15
css_cluster_distance_method <- "pearson" ## suggested to make comparable across samples
# css_cluster_distance_method <- "euclidean" ## not suggested
css_cluster_results <- ConsensusClusterPlus(
    d = mat_score,
    maxK = maxK, # Maximum number of clusters to evaluate
    reps = 100, # Number of resampling iterations
    pItem = 0.8, # Proportion of items to sample
    pFeature = 1, # Proportion of features to sample
    clusterAlg = "hc", # Clustering algorithm (hierarchical clustering)
    distance = css_cluster_distance_method, # or euclidean or pearson Distance metric
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
if (file.exists(file.path(dir_CSSCluster, "CDF_values.rds"))) {
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
if (F) {
    pdf(file.path(dir_CSSCluster, "NbClust_GapStat.pdf"), width = 4, height = 4)
    ## takes 1min 
    tmp <- mat_score
    tmp[is.na(tmp)] <- 0
    p <- fviz_nbclust(x = t(tmp), FUNcluster = hcut, method = "gap_stat", k.max = maxK)
    rm(tmp)
    print(p)
    dev.off()
}
pdf(file.path(dir_CSSCluster, "NbClust_silhouette.pdf"), width = 4, height = 4)
p <- fviz_nbclust(x = t(mat_score), FUNcluster = hcut, method = "silhouette", k.max = maxK)
print(p)
dev.off()

#------------------ ~~~ Heatmap with the CSS clusters ~~~ --------------------
cli_h1("Heatmap with the CSS clusters")

df_css_cluster <- lapply(2:maxK, function(k) {
    o <- data.frame(
        ID = names(css_cluster_results[[k]]$consensusClass),
        res = css_cluster_results[[k]]$consensusClass
    )
    colnames(o)[2] <- paste0("cut", k)
    o
}) %>% Reduce(left_join, .)
# view(df_css_cluster)
column_to_rownames(df_css_cluster, "ID") -> df_css_cluster
write_rds(df_css_cluster, file.path(dir_res, glue("dataframe.cellmeta_css_clusters_{res_str}.rds")))

pal_css_cluster <- lapply(2:maxK, function(k) {
    cl <- css_cluster_results[[k]]$clrs[[3]]
    if (length(cl) < k) {
        cl <- rainbow(n = k)
    }
    o <- structure(cl, names = 1:length(cl))
    o
})
names(pal_css_cluster) <- paste0("cut", 2:maxK)
write_rds(pal_css_cluster, file.path(dir_res, glue("palette_css_clusters_{res_str}.rds")))
# pal_css_cluster <- readRDS(file.path(dir_res, glue("palette_css_clusters_{res_str}.rds")))

ht_col_order <- css_cluster_results[[maxK]]$consensusTree$order

#------------------ ~~~ Try a cut ~~~ --------------------
cli_h1("Try a cut")
tentative_cut <- 8
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
    mat_score[, ht_col_order],
    name = "zscore",
    top_annotation = ha_top,
    cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
    show_column_names = F, cluster_columns = TRUE, show_column_dend = FALSE, column_names_side = "top", column_dend_height = unit(5, "mm"),
    column_split = tentative_cut, column_gap = unit(0, "mm"), border = TRUE,
    column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
    column_title = res_str,
    col = heatmap_color_fun_zscore_2, na_col = "ghostwhite",
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
    width = unit(7, "inches"),
    height = unit(5, "inches"), raster_by_magick = TRUE
)
pdf(file.path(dir_res, glue("heatmap_niche_score_matrix_all_samples_css_clusters_{res_str}.hcut.pdf")),
    width = 7 + 1, height = 5 + 1 + 0.08 * maxK, useDingbats = FALSE
)
draw(ht_css, heatmap_legend_side = "bottom")
dev.off()

ht_css <- ComplexHeatmap::Heatmap(
    mat_score[, ht_col_order],
    name = "zscore",
    top_annotation = ha_top,
    cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
    show_column_names = F, cluster_columns = FALSE, show_column_dend = FALSE, column_names_side = "top", column_dend_height = unit(5, "mm"),
    # column_split = tentative_cut, column_gap = unit(0, "mm"), border = TRUE,
    column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
    column_title = res_str,
    col = heatmap_color_fun_zscore_2, na_col = "ghostwhite",
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
    width = unit(7, "inches"),
    height = unit(5, "inches"), raster_by_magick = TRUE
)
pdf(file.path(dir_res, glue("heatmap_niche_score_matrix_all_samples_css_clusters_{res_str}.uncut.pdf")),
    width = 7 + 1, height = 5 + 1 + 0.08 * maxK, useDingbats = FALSE
)
draw(ht_css, heatmap_legend_side = "bottom")
dev.off()

#------------------ ~~~ Propose the optimal K ~~~ --------------------
cli_h1("Propose the optimal K")
for (optimal_cut in c(8, 9, 10, 11, 12)) {
    cli_rule(glue("Try cut={optimal_cut}"))
    f_dict_niche_to_metaniche <- file.path(dir_res, glue("dict_niche_to_metaniche_cut{optimal_cut}.csv"))

    dict_niche_to_metaniche <- structure(
        df_css_cluster[[paste0("cut", optimal_cut)]],
        names = rownames(df_css_cluster)
    )
    dict_niche_to_metaniche <- df_css_cluster[, paste0("cut", optimal_cut), drop = F] %>%
        rownames_to_column("niche")
    colnames(dict_niche_to_metaniche) <- c("niche", "MetaNiche")
    print(head(dict_niche_to_metaniche))
    dict_niche_to_metaniche$MetaNiche <- paste0("N", dict_niche_to_metaniche$MetaNiche)
    metaniche_str <- "MetaNiche"
    write.csv(dict_niche_to_metaniche, f_dict_niche_to_metaniche)
    write_rds(dict_niche_to_metaniche, file.path(dir_res, glue("dict_niche_to_metaniche_cut{optimal_cut}.rds")))
    if (T) {
        ## Re-plot the heatmap with the optimal cut
        ht_col_order <- css_cluster_results[[maxK]]$consensusTree$order ## ensure the order is the same
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
            mat_score[, ht_col_order],
            name = "zscore",
            top_annotation = ha_top,
            cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
            row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
            show_column_names = F, cluster_columns = TRUE, cluster_column_slices = FALSE,
            show_column_dend = FALSE, column_names_side = "top", column_dend_height = unit(5, "mm"),
            column_split = df_css_cluster[ht_col_order, paste0("cut", optimal_cut)], 
            column_gap = unit(0, "mm"), border = TRUE,
            column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
            column_title = res_str,
            col = heatmap_color_fun_zscore_2, na_col = "ghostwhite",
            heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
            width = unit(7, "inches"),
            height = unit(5, "inches"), raster_by_magick = TRUE
        )
        pdf(file.path(dir_res, glue("heatmap_niche_score_matrix_all_samples_css_clusters_{res_str}.cut{optimal_cut}.pdf")),
            width = 7 + 1, height = 5 + 1 + 0.08 * maxK, useDingbats = FALSE
        )
        draw(ht_css, heatmap_legend_side = "bottom")
        dev.off()
    }


    # f_metaniche_assay <- file.path(dir_res, glue("metaniche_assay_css_clusters_{res_str}_cut{optimal_cut}.rds"))
    f_metaniche_assay <- file.path(dir_res, glue("metaniche_assay_cut{optimal_cut}.rds"))
    if (file_exists(f_metaniche_assay)) {
        metaniche_assay <- readRDS(f_metaniche_assay)
    } else {
        #------------------ ~~~ Create the MetaNiche assay ~~~ --------------------
        cli_h1("Create the MetaNiche assay")

        mat_metaniche_relative <- mat_tapply(t(mat_data),
            INDEX = dict_niche_to_metaniche$MetaNiche,
            FUN = mean
        ) %>% t()
        # sum(mat_metaniche[, 1]) == 100
        metaniche_assay <- CreateAssayObject(
            data = mat_metaniche_relative, assay = "niche", min.cells = 0,
            min.features = 0
        )
        metaniche_assay <- ScaleData(
            metaniche_assay,
            features = rownames(metaniche_assay),
            do.center = TRUE, do.scale = TRUE, scale.max = 4
        )
        metaniche_assay <- CreateSeuratObject(
            metaniche_assay,
            assay = "niche", meta.data = NULL
        )
        write_rds(metaniche_assay, f_metaniche_assay)
    }
    print(metaniche_assay)

    write_rds(
        GetAssayData(metaniche_assay, assay = "niche", layer = "data"),
        file.path(dir_res, glue("metaniche_data_matrix_cut{optimal_cut}.rds"))
    )
    write_rds(
        GetAssayData(metaniche_assay, assay = "niche", layer = "scale.data"),
        file.path(dir_res, glue("metaniche_scale_matrix_cut{optimal_cut}.rds"))
    )

    pal_metaniche <- pal_css_cluster[[paste0("cut", optimal_cut)]]
    names(pal_metaniche) <- paste0("N", names(pal_metaniche))
    write_rds(pal_metaniche, file.path(dir_res, glue("palette_metaniche_cut{optimal_cut}.rds")))

    #------------------ ~~~ Heatmap the metaniche ~~~ --------------------
    cli_h1("Heatmap the metaniche")
    mat_viz <- GetAssayData(metaniche_assay, assay = "niche", layer = "scale.data")
    print(range(mat_viz))
    ht_col_anno <- ComplexHeatmap::HeatmapAnnotation(
        MetaNiche = colnames(mat_viz),
        col = list(MetaNiche = pal_metaniche),
        which = "column",
        show_legend = F, show_annotation_name = F, annotation_name_gp = gpar(fontsize = 6),
        simple_anno_size = unit(2, "mm"), gap = 0
    )
    ht_metaniche <- ComplexHeatmap::Heatmap(
        mat_viz,
        name = "zscore",
        top_annotation = ht_col_anno,
        cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
        row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
        show_column_names = TRUE, cluster_columns = TRUE, show_column_dend = TRUE, column_names_side = "top", column_dend_height = unit(3, "mm"),
        column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
        column_title = metaniche_str,
        col = heatmap_color_fun_zscore_2, na_col = "ghostwhite",
        heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
        width = unit(5 * ncol(mat_viz) / nrow(mat_viz), "inches"),
        height = unit(5, "inches"), raster_by_magick = TRUE
    )
    pdf(file.path(dir_res, glue("heatmap_metaniche_scale_matrix_cut{optimal_cut}.pdf")),
        width = 5 * ncol(mat_viz) / nrow(mat_viz) + 1,
        height = 5 + 2, useDingbats = FALSE
    )
    draw(ht_metaniche, heatmap_legend_side = "bottom")
    dev.off()

    mat_viz <- GetAssayData(metaniche_assay, assay = "niche", layer = "data")
    ht_col_anno <- ComplexHeatmap::HeatmapAnnotation(
        MetaNiche = colnames(mat_viz),
        col = list(MetaNiche = pal_metaniche),
        which = "column",
        show_legend = F, show_annotation_name = F, annotation_name_gp = gpar(fontsize = 6),
        simple_anno_size = unit(2, "mm"), gap = 0
    )
    ht_metaniche <- ComplexHeatmap::Heatmap(
        mat_viz,
        name = "zscore",
        top_annotation = ht_col_anno,
        cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left",
        row_names_gp = gpar(fontsize = 6), clustering_method_rows = "ward.D2",
        show_column_names = TRUE, cluster_columns = TRUE, show_column_dend = TRUE, column_names_side = "top", column_dend_height = unit(3, "mm"),
        column_names_gp = gpar(fontsize = 6), clustering_method_columns = "ward.D2",
        column_title = metaniche_str,
        na_col = "ghostwhite",
        heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
        width = unit(5 * ncol(mat_viz) / nrow(mat_viz), "inches"),
        height = unit(5, "inches"), raster_by_magick = TRUE
    )
    pdf(file.path(dir_res, glue("heatmap_metaniche_data_matrix_cut{optimal_cut}.pdf")),
        width = 5 * ncol(mat_viz) / nrow(mat_viz) + 1,
        height = 5 + 2, useDingbats = FALSE
    )
    draw(ht_metaniche, heatmap_legend_side = "bottom")
    dev.off()
}

#------------------ ~~~ Check if metaniche are patient specific ~~~ --------------------
for (optimal_cut in c(8, 9, 10, 11, 12)) {
    cli_rule(glue("Try cut={optimal_cut}"))
    dict_niche_to_metaniche <- read_rds(file.path(dir_res, glue("dict_niche_to_metaniche_cut{optimal_cut}.rds")))
    cli_h1("Check if metaniche are patient specific")
    dict_niche_to_metaniche <- dict_niche_to_metaniche %>% 
        separate(niche, into = c("Sample", "Niche"), sep = "_", remove = FALSE)
    metaniche_n_samples <- dict_niche_to_metaniche %>%
        group_by(MetaNiche) %>%
        summarise(n_sample = n_distinct(Sample)) %>%
        arrange(desc(MetaNiche))
    metaniche_n_samples$MetaNiche <- factor(
        metaniche_n_samples$MetaNiche, 
        levels = gtools::mixedsort(unique(metaniche_n_samples$MetaNiche)))
    ## barplot n_samples across metaniches
    p <- ggbarplot(
        metaniche_n_samples, x = "MetaNiche", y = "n_sample",
        fill = "MetaNiche", palette = pal_metaniche,
        label = TRUE, lab.size = 3, lab.pos = "out",
        xlab = "MetaNiche", 
        ylab = sprintf("Number of samples (N total = %s)", length(sample_array)),
        # sort.val = "desc", sort.by.groups = FALSE,
        legend = "none"
    )
    ggsave(
        filename = file.path(dir_res, glue("barplot_metaniche_n_samples_cut{optimal_cut}.pdf")),
        plot = p, width = 4, height = 3)
}

#------------------ ~~~ Attach the metaniche info to each samples ~~~ --------------------
cli_h1("Attach the metaniche info to each samples")
for (optimal_cut in c(10)) {
    pal_metaniche <- read_rds(file.path(dir_res, glue("palette_metaniche_cut{optimal_cut}.rds")))
    cli_rule(glue("Try cut={optimal_cut}"))
    dict_niche_to_metaniche <- read_rds(file.path(dir_res, glue("dict_niche_to_metaniche_cut{optimal_cut}.rds")))
    dict_niche_to_metaniche <- dict_niche_to_metaniche %>% 
        separate(niche, into = c("Sample", "Niche"), sep = "_", remove = FALSE)

    for (s in sample_array) {
        cli_h2(s)
        dir_compo_s <- file.path(dir_proj, "split_sample", s, "niche_proposal_onestep/snn_clusters")
        dict_niche_to_metaniche_s <- dict_niche_to_metaniche %>% filter(Sample == s)

        f_cellmeta <- file.path(dir_compo_s, "dataframe.cellmeta.allcolumns.with_niche_cluster.rds")
        f_cellmeta_out <- file.path(dir_compo_s, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}.rds"))

        # if (!file_exists(f_cellmeta)) {
        #     cli_alert_danger("No cell meta data file!")
        #     next()
        # }
        # if (file_exists(f_cellmeta_out)) {
        #     cli_alert_info("Output file already exists, skip.")
        #     next()
        # }

        if (file_exists(f_cellmeta) & !file_exists(f_cellmeta_out)) {
        # if (T) {
            df_cellmeta_s <- readRDS(f_cellmeta)
            # print(head(df_cellmeta_s))
            df_cellmeta_s$Niche <- paste0('g', df_cellmeta_s$niche_cluster)
            df_cellmeta_s <- left_join(df_cellmeta_s, dict_niche_to_metaniche_s[, c("Niche", "MetaNiche")], by = "Niche")
            WANTED_COLNAMES <- c(
                "celltype", "cellname", "sample", 
                "coord_x", "coord_y", 
                "coord_x_adj_to_origin", "coord_y_adj_to_origin", 
                "coord_x_adj", "coord_y_adj", 
                "cell_state_paper", "pCR_status",
                "archetype", "MetaNiche"
            )
            df_cellmeta_s <- df_cellmeta_s[, intersect(WANTED_COLNAMES, colnames(df_cellmeta_s))]
            rownames(df_cellmeta_s) <- df_cellmeta_s$cellname
            write_rds(df_cellmeta_s, f_cellmeta_out)
            write_csv(df_cellmeta_s, file.path(dir_compo_s, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}.csv")))
            print(head(df_cellmeta_s))
        }
    }
    cat("[done] for cut=", optimal_cut, "\n", sep = "")
}

#------------------ ~~~ Conver the metaniche data frame to Xenium browser ~~~ --------------------
cli_h1("Conver the metaniche data frame to Xenium browser")
for (optimal_cut in c(10)) {
    for (s in sample_array) {
        cli_h2(s)
        dir_compo_s <- file.path(dir_proj, "split_sample", s, "niche_proposal_onestep/snn_clusters")
        f_cellmeta <- file.path(dir_compo_s, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}.rds"))
        f_xenium_browser <- file.path(dir_compo_s, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}_for_xenium_browser.csv"))
        if (!file_exists(f_xenium_browser)) {
            df_cellmeta_s <- readRDS(f_cellmeta)
            df_xenium_browser <- df_to_exnium_explorer(df_cellmeta_s, 'cellname', 'MetaNiche')
            write_csv(df_xenium_browser, f_xenium_browser)
        }
    }
}

#------------------ ~~~ Combine the single-cell metaniche data frames ~~~ --------------------
cli_h1("Combine the single-cell metaniche data frames")
for (optimal_cut in c(10)) {
    pal_metaniche <- read_rds(file.path(dir_res, glue("palette_metaniche_cut{optimal_cut}.rds")))
    cli_rule(glue("Merge cut={optimal_cut}"))
    list_df_cellmeta_metaniche <- lapply(sample_array, function(s) {
        dir_compo_s <- file.path(dir_proj, "split_sample", s, "niche_proposal_onestep/snn_clusters")
        f_cellmeta_metaniche <- file.path(dir_compo_s, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}.rds"))
        if (!file_exists(f_cellmeta_metaniche)) {
            cli_alert_danger("No cell meta data file!")
            return(NULL)
        }
        df_cellmeta_s <- readRDS(f_cellmeta_metaniche)
        rownames(df_cellmeta_s) <- paste0(s, "_", df_cellmeta_s$cellname)
        return(df_cellmeta_s)
    })
    names(list_df_cellmeta_metaniche) <- sample_array
    df_cellmeta_metaniche_all <- bind_rows(list_df_cellmeta_metaniche)

    df_cellmeta_metaniche_all$MetaNiche <- factor(
        df_cellmeta_metaniche_all$MetaNiche,
        levels = gtools::mixedsort(unique(df_cellmeta_metaniche_all$MetaNiche))
    )
    write_rds(
        df_cellmeta_metaniche_all,
        file.path(dir_res, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}_all_samples.rds"))
    )
    write_csv(
        df_cellmeta_metaniche_all,
        file.path(dir_res, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}_all_samples.csv"))
    )
    print(head(df_cellmeta_metaniche_all, 3))
    print(tail(df_cellmeta_metaniche_all, 3))
    print(table(df_cellmeta_metaniche_all$MetaNiche))
}

for (optimal_cut in c(10)) {
    pal_metaniche <- read_rds(file.path(dir_res, glue("palette_metaniche_cut{optimal_cut}.rds")))
    ggsave(file.path(dir_res, glue("palette_metaniche_cut{optimal_cut}.pdf")),
        pal_to_ggplot(pal_metaniche, "MetaNiche"), 
        width = 2, height = length(pal_metaniche)*0.5, useDingbats = FALSE)
}

#------------------ ~~~ Spatial visualize meta niche for each sample ~~~ --------------------
cli_h1("Spatial visualize meta niche for each sample")
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
fs::dir_create(file.path(dir_res, "spatial_plots_each_sample"))
for (optimal_cut in c(10)) {
    pal_metaniche <- read_rds(file.path(dir_res, glue("palette_metaniche_cut{optimal_cut}.rds")))
    cli_rule(glue("Spatial Viz cut={optimal_cut}"))
    for (s in sample_array) {
        cli_h2(s)
        dir_compo_s <- file.path(dir_proj, "split_sample", s, "niche_proposal_onestep/snn_clusters")
        f_cells_metaniche_s <- file.path(dir_compo_s, glue("dataframe.cellmeta.allcolumns.with_metaniche_cut{optimal_cut}.rds"))
        df_cellmeta_s <- readRDS(f_cells_metaniche_s)
        
        f_xenium <- glue("/volumes/USR1/yyan/project/tnbc_xenium/data/{s}/cleaned1/ready.seurat.rds")
        xmo <- readRDS(f_xenium)
        # head(df_cellmeta_s)
        all(rownames(df_cellmeta_s) %in% colnames(xmo))
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

        #------ SpatialDimPlot ------
        pal_metaniche["NA"] <- "ghostwhite"
        metaniche_str <- "MetaNiche"
        xmo <- AddMetaData(xmo, metadata = df_cellmeta_s[colnames(xmo), metaniche_str, drop = TRUE], col.name = metaniche_str)
        p <- ImageDimPlotSegmentation.raster(
            xmo,
            group.by = metaniche_str, cols = pal_metaniche,
            cell_segmentation_width = 0.01
        ) +
            labs(title = glue("Meta Niche(k={optimal_cut})")) +
            rremove("legend")
        p <- q_add_anno_xenium_scale_bar(
            p,
            roi_xmin = tissue_xmin, roi_xmax = tissue_xmax, roi_ymin = tissue_ymin, roi_ymax = tissue_ymax,
            color = ifelse(insitu_draw_dark_bg, "white", "black"), linewidth = 3
        )
        ggsave(
            filename = file.path(
                dir_res, "spatial_plots_each_sample",
                glue("spatial_dimplot_{metaniche_str}_{optimal_cut}_{s}.pdf")
            ),
            plot = p,
            width = tissue_img_width, height = tissue_img_height, useDingbats = FALSE
        )
        cat("[done] for sample=", s, "\n", sep = "")
    }
}
#------------------ ~~~ Export ~~~ --------------------
cli_h1("Export")
cat("[done]")
timestamp()
