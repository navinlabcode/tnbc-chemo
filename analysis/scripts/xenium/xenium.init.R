#!/usr/bin/Rscript
suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal:init xenium seurat object
# R: 4.4.1
# Seurat: 5.1.0
# Known error:  File '/volumes/USR1/yyan/project/tnbc_xenium/data0/ART94/transcripts.csv.gz' does not exist or is non-readable.
# transcripts.parquet is what it saved
# 
# Possible solution: Install https://github.com/10XGenomics/seurat/tree/develop or use my previous script xenium.init.Seurat4.R.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-08-17 / 2024-12-23
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
  library(ruok)
  library(Seurat)
  # library(arrow)
  library(fs)
  library(Signac)
  library(clustree)
})
options <- commandArgs(trailingOnly = TRUE)

#------------------ ~~~ input ~~~ --------------------
sample_name <- "ART10"
do_binarize_counts <- F
do_rna_like_preprocessing <- T
do_atac_like_preprocessing <- F

if (length(options) > 0) {
  sample_name <- options[[1]]
  do_binarize_counts <- as.numeric(options[[2]]) == 1
  do_rna_like_preprocessing <- as.numeric(options[[3]]) == 1
  do_atac_like_preprocessing <- as.numeric(options[[4]]) == 1
}

n_linear_dr_components <- 50
n_linear_dr_components_for_umap <- 50

#------ find input and setup output ------
dir_in <- file.path(
  "/volumes/USR1/yyan/project/tnbc_xenium",
  "data0", sample_name
)
cli_alert_info(c("Processing ", sample_name, " ", dir_in))
# stopifnot(dir_exists(dir_in))

dir_res <- file.path(
  "/volumes/USR1/yyan/project/tnbc_xenium",
  "data",
  sample_name
)
dir_create(dir_res)

if (sum(c(do_rna_like_preprocessing, do_atac_like_preprocessing)) != 1) {
  stop("Either RNA-like or ATAC-like preprocessing is allowed. Cannot do both or none. ")
}
#------------------ ~~~ helper funcs ~~~ --------------------

# # Redefine ReadXenium()
# ReadXenium.adhoc <- function(
#     data.dir, outs = c("matrix", "microns"),
#     type = "centroids",
#     mols.qv.threshold = 20) {
#   type <- match.arg(
#     arg = type, choices = c("centroids", "segmentations"),
#     several.ok = TRUE
#   )
#   outs <- match.arg(
#     arg = outs, choices = c("matrix", "microns"),
#     several.ok = TRUE
#   )
#   outs <- c(outs, type)
#   has_dt <- requireNamespace("data.table", quietly = TRUE) &&
#     requireNamespace("R.utils", quietly = TRUE)
#   data <- sapply(outs, function(otype) {
#     switch(EXPR = otype,
#       matrix = {
#         matrix <- suppressWarnings(Read10X(data.dir = file.path(
#           data.dir,
#           "cell_feature_matrix/"
#         )))
#         matrix
#       },
#       centroids = {
#         if (has_dt) {
#           cell_info <- as.data.frame(data.table::fread(file.path(
#             data.dir,
#             "cells.csv.gz"
#           )))
#         } else {
#           cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
#         }
#         cell_centroid_df <- data.frame(
#           x = cell_info$x_centroid,
#           y = cell_info$y_centroid, cell = cell_info$cell_id,
#           stringsAsFactors = FALSE
#         )
#         cell_centroid_df
#       },
#       segmentations = {
#         if (has_dt) {
#           cell_boundaries_df <- as.data.frame(data.table::fread(file.path(
#             data.dir,
#             "cell_boundaries.csv.gz"
#           )))
#         } else {
#           cell_boundaries_df <- read.csv(file.path(
#             data.dir,
#             "cell_boundaries.csv.gz"
#           ), stringsAsFactors = FALSE)
#         }
#         names(cell_boundaries_df) <- c("cell", "x", "y")
#         cell_boundaries_df
#       },
#       microns = {
#         transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
#         transcripts <- subset(transcripts, qv >= mols.qv.threshold)

#         df <- data.frame(
#           x = transcripts$x_location, y = transcripts$y_location,
#           gene = transcripts$feature_name, stringsAsFactors = FALSE
#         )
#         df
#       },
#       stop("Unknown Xenium input type: ", otype)
#     )
#   }, USE.NAMES = TRUE)
#   return(data)
# }


#------------------ ~~~ Rock ~~~ --------------------

#------ import data ------

# if (!file_exists(file.path(dir_res, "xenium.seurat.rds"))) {
#   data <- ReadXenium.adhoc(dir_in,
#     outs = c("matrix", "microns"),
#     type = c("centroids", "segmentations")
#   )
#   names(data)
#   ## continue the regular LoadXenium
#   segmentations.data <- list(
#     centroids = CreateCentroids(data$centroids),
#     segmentation = CreateSegmentation(data$segmentations)
#   )
#   coords <- CreateFOV(
#     coords = segmentations.data,
#     type = c("segmentation", "centroids"),
#     molecules = data$microns,
#     assay = "Xenium"
#   )
#   xmo <- CreateSeuratObject(
#     counts = data$matrix[["Gene Expression"]],
#     assay = "Xenium"
#   )
#   xmo[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
#   xmo[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
#   xmo[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
#   xmo[["fov"]] <- coords

# if (!file_exists(file.path(dir_res, "xenium.seurat.rds"))) {
if (T) {
  xmo <- LoadXenium(dir_in, fov = "fov", assay = "Xenium")
  print(xmo)

  # An object of class Seurat
  # 6717 features across 39525 samples within 4 assays
  # Active assay: Xenium (5001 features, 0 variable features)
  # 3 other assays present: BlankCodeword, ControlCodeword, ControlProbe
  # 1 spatial field of view present: fov

  # str(rownames(xmo))
  # chr [1:5001] "A2ML1" "AAMP" "AAR2" "AARSD1" "ABAT" "ABCA1" "ABCA3"

  write_rds(xmo, file.path(dir_res, "xenium0.seurat.rds"))

  #------ remove cells with 0 counts ------
  n_cells_raw <- ncol(xmo)
  n_cells_nonempty <- sum(xmo$nCount_Xenium != 0)

  write_lines(
    x = sprintf(
      "%s / %s (%.3f%%) cells have non-zero counts",
      comma(n_cells_nonempty),
      comma(n_cells_raw),
      n_cells_nonempty / n_cells_raw * 100
    ),
    file = file.path(dir_res, "log.qc.remove_cells_progress.txt")
  )

  xmo <- subset(xmo, cells = Cells(xmo)[xmo$nCount_Xenium != 0])
  print(xmo)
  write_rds(xmo, file.path(dir_res, "xenium.seurat.rds"))
} else {
  xmo <- read_rds(file.path(dir_res, "xenium.seurat.rds"))
}

#------ depending on preprocessing procedures ------
if (do_binarize_counts & do_rna_like_preprocessing) {
  dir_res <- file.path(dir_res, "binarized_pca")
  f_out <- file.path(dir_res, "xenium_binarized_pca.seurat.rds")
}
if (do_binarize_counts & do_atac_like_preprocessing) {
  # most common for scATAC
  dir_res <- file.path(dir_res, "binarized_svd")
  f_out <- file.path(dir_res, "xenium_binarized_svd.seurat.rds")
}
if (!do_binarize_counts & do_rna_like_preprocessing) {
  # most common for scRNA
  dir_res <- file.path(dir_res, "nonbinarized_pca")
  f_out <- file.path(dir_res, "xenium_nonbinarized_pca.seurat.rds")
}
if (!do_binarize_counts & do_atac_like_preprocessing) {
  # not recommended
  dir_res <- file.path(dir_res, "nonbinarized_svd")
  f_out <- file.path(dir_res, "xenium_nonbinarized_svd.seurat.rds")
}
dir_create(dir_res)
message(dir_res)
message(f_out)

# if (file_exists(f_out)) {
if (F) {
  message("Load the existing object...")
  xmo <- read_rds(f_out)
} else {
  message("Processing begins...")
  #------------------ ~~~ Preprocessing ~~~ --------------------
  if (do_binarize_counts) {
    cli_alert_warning("Binarizing count as requested...")
    xmo <- BinarizeCounts(xmo, assay = "Xenium")
    ## --- to improve ---
    ## Bug: BinarizeCounts will re-create the nCount and nFeature based on the binarized result.
  }
  write_rds(xmo, f_out)

  if (sum(c(do_rna_like_preprocessing, do_atac_like_preprocessing)) != 1) {
    stop("Either RNA-like or ATAC-like preprocessing is allowed. Cannot do both or none. ")
  }
  if (do_rna_like_preprocessing) {
    cli_alert_info("Running RNA-like processing with PCA")
  }
  if (do_atac_like_preprocessing) {
    cli_alert_info("Running ATAC-like processing with LDA")
  }

  name_linear_dr <- NULL
  if (do_rna_like_preprocessing) {
    name_linear_dr <- "pca"
  }
  if (do_atac_like_preprocessing) {
    name_linear_dr <- "lsi"
  }

  i_linear_dr_components_for_umap <- 1:n_linear_dr_components_for_umap

  DefaultAssay(xmo)
  #------ linear dimension reduction ------
  ## PCA or TFIDF
  if (do_atac_like_preprocessing) {
    xmo <- RunTFIDF(xmo, assay = "Xenium")
    xmo <- FindTopFeatures(
      xmo,
      min.cutoff = 0, assay = "Xenium"
    )
    str(VariableFeatures(xmo))
    xmo <- RunSVD(xmo, n = n_linear_dr_components)

    p <- DepthCor(xmo, assay = "Xenium", reduction = name_linear_dr, n = n_linear_dr_components_for_umap) +
      geom_hline(yintercept = c(-0.75, 0.75), lty = "dashed", color = "red")
    ruok::ggsave2(file.path(dir_res, "qc.DepthCor"), p, width = 5, height = 3)

    i_linear_dr_components_for_umap <- setdiff(
      i_linear_dr_components_for_umap,
      which(abs(p$data$counts) > 0.75)
    )
    cat("these LSI components are used: \n")
    print(i_linear_dr_components_for_umap)
    print(xmo)
  }
  if (do_rna_like_preprocessing) {
    is_non_empty_features <- rowSums(xmo@assays$Xenium@counts) != 0
    table(is_non_empty_features)
    VariableFeatures(xmo) <- rownames(xmo@assays$Xenium@counts)[is_non_empty_features]
    xmo <- NormalizeData(xmo)
    xmo <- ScaleData(xmo)
    xmo <- RunPCA(xmo, npcs = n_linear_dr_components)
  }

  write_rds(xmo, f_out)

  #------ umap ------
  xmo <- RunUMAP(xmo, reduction = name_linear_dr, dims = i_linear_dr_components_for_umap)
  write_rds(xmo, f_out)
  #
  # p <- FeaturePlot(
  #   xmo,
  #   features = 'nCount_Xenium',
  #   reduction = 'umap',
  #   max.cutoff = 'q95',
  #   min.cutoff = 'q5'
  # )
  # p
  # ruok::ggsave2(file.path(dir_res, 'umap.nCount_Xenium'), p, width=7, height = 7)

  #------ cell clustering ------
  xmo <- FindNeighbors(
    xmo,
    k.param = 20,
    dims = i_linear_dr_components_for_umap,
    reduction = name_linear_dr
  )
  write_rds(xmo, f_out)


  idx <- grepl(pattern = "_snn_res", x = colnames(xmo[[]]))
  for (i in colnames(xmo[[]])[idx]) {
    xmo[[i]] <- NULL
  }
  snn_res_max <- 0.8
  if (ncol(xmo) < 100) {
    snn_res_max <- 1
  }

  for (snn_res_i in seq(from = 0.2, to = snn_res_max, by = .2)) {
    cat(snn_res_i, "... ")
    try(xmo <- FindClusters(xmo,
      resolution = snn_res_i,
      algorithm = 3,
      verbose = F
    ))
    snn_str_i <- sprintf("%s_snn_res.%s", DefaultAssay(xmo), snn_res_i)
    if (snn_str_i %in% colnames(xmo[[]])) {
      xmo[[snn_str_i]] <- Idents(xmo) ## The cluster order is human-readable
    }
  }
  cat("\n")
  write_rds(xmo, f_out)

  p <- clustree::clustree(
    xmo,
    prefix = sprintf("%s_snn_res.", DefaultAssay(xmo)),
    prop_filter = 0, layout = "sugiyama"
  )
  ruok::ggsave2(file.path(dir_res, sprintf("qc.clustree")),
    p,
    width = 15, height = 9
  )
}





#------------------ ~~~ Visualization ~~~ --------------------
theme_set(theme_pubr(base_size = 7, legend = "right"))
cli_h1("Visualization")
fov <- xmo[["fov"]]
df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
df_coord <- df_coord[Cells(xmo), ]
coord_xy_ratio <- diff(range(df_coord$x)) / diff(range(df_coord$y))
if (coord_xy_ratio < 1) {
  pdf_img_width <- 7.5
  pdf_img_height <- pdf_img_width * coord_xy_ratio
  pdf_img_height <- pmax(pdf_img_height, 3)
} else {
  pdf_img_height <- 7.5
  pdf_img_width <- pdf_img_height / coord_xy_ratio
  pdf_img_width <- pmax(pdf_img_width, 3)
}

#------ viz technical features ------

for (feature_i in c("nCount_Xenium", "nFeature_Xenium", "EPCAM")) {
  message(feature_i)
  p <- ImageFeaturePlot(
    xmo,
    features = feature_i,
    size = 0.5,
    max.cutoff = "q95",
    min.cutoff = "q5", dark.background = F,
    cols = c("lightgrey", "maroon1")
  )
  ruok::ggsave2(
    file.path(
      dir_res,
      sprintf("imgfeature.%s", feature_i)
    ),
    p,
    width = pdf_img_width, height = pdf_img_height
  )
  p <- FeaturePlot(
    xmo,
    features = feature_i,
    reduction = "umap",
    max.cutoff = "q95",
    min.cutoff = "q5", order = T, raster = T, raster.dpi = c(1024, 1024),
    cols = c("lightgrey", "maroon1")
  )
  ruok::ggsave2(
    file.path(
      dir_res,
      sprintf("featureplot.%s", feature_i)
    ),
    p,
    width = 7, height = 7
  )
}

#------ viz idents ------

snn_res_i <- 0.2
snn_str_i <- sprintf("%s_snn_res.%s", DefaultAssay(xmo), snn_res_i)

snn_str_opts <- colnames(xmo@meta.data)[grepl(sprintf("%s_snn_res", DefaultAssay(xmo)), colnames(xmo@meta.data))]
for (snn_str_i in snn_str_opts) {
  xmo$seurat_clusters <- xmo@meta.data[, snn_str_i]

  pal_snn <- structure(
    Seurat::DiscretePalette(n = length(levels(xmo$seurat_clusters)), palette = "parade"),
    names = levels(xmo$seurat_clusters)
  )
  pa <- DimPlot(
    xmo,
    reduction = "umap", label = T,
    group.by = snn_str_i, shuffle = T,
    cols = pal_snn
  )

  ruok::ggsave2(
    file.path(
      dir_res,
      sprintf("dimplot.%s", snn_str_i)
    ),
    pa,
    width = 7, height = 7
  )

  pb <- ImageDimPlot(
    object = xmo, cols = "parade",
    group.by = "seurat_clusters", dark.background = FALSE
  )

  ruok::ggsave2(
    file.path(
      dir_res,
      sprintf("imgdimplot.%s", snn_str_i)
    ),
    pb,
    width = pdf_img_width, height = pdf_img_height
  )
}

#------ animation of UMAP to spacial coordinates ------

library(scattermore)
library(gganimate)

for (snn_str_i in snn_str_opts) {
  xmo$seurat_clusters <- xmo@meta.data[, snn_str_i]

  pal_snn <- structure(
    Seurat::DiscretePalette(n = length(levels(xmo$seurat_clusters)), palette = "parade"),
    names = levels(xmo$seurat_clusters)
  )

  fov <- xmo[["fov"]]
  head(fov$centroids@coords)
  head(fov$centroids@cells)
  df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
  df_coord <- df_coord[Cells(xmo), ]
  coord_xy_ratio <- diff(range(df_coord$x)) / diff(range(df_coord$y))
  df_coord$x <- rescale(df_coord$x, to = c(-1 * coord_xy_ratio, coord_xy_ratio))
  df_coord$y <- -1 * rescale(df_coord$y, to = c(-1, 1))
  df_spatial <- data.frame(df_coord, seurat_clusters = xmo$seurat_clusters, type = "spatial", stringsAsFactors = F)

  df_umap <- as.data.frame(Embeddings(xmo, "umap"))
  colnames(df_umap) <- c("x", "y")
  df_umap$x <- rescale(df_umap$x, to = c(-1, 1))
  df_umap$y <- rescale(df_umap$y, to = c(-1, 1))
  df_umap <- data.frame(df_umap, seurat_clusters = xmo$seurat_clusters, type = "embedding", stringsAsFactors = F)

  tmp <- sample(1:nrow(df_spatial), size = round(nrow(df_spatial) / 3), replace = F)
  df_trans <- rbind(df_spatial[tmp, ], df_umap[tmp, ])
  # df_trans <- rbind(df_spatial, df_umap)

  panim <- ggplot(df_trans, aes(x = x, y = y)) +
    geom_scattermore(
      aes(x = x, y = y, color = seurat_clusters),
      pointsize = .5
    ) +
    scale_color_manual(values = pal_snn) +
    Seurat::DarkTheme() +
    rremove("x.title") +
    rremove("y.title") +
    rremove("x.text") +
    rremove("y.text") +
    rremove("x.axis") +
    rremove("y.axis") +
    coord_equal() +
    theme(legend.position = "none")
  # panim

  panim <- panim +
    transition_states(type,
      transition_length = 5,
      state_length = 1
    ) +
    labs(title = "{closest_state}") +
    theme(plot.title = element_text(size = 28)) +
    enter_fade()

  anim_res <- panim + view_follow()
  anim_save(file.path(
    dir_res,
    sprintf("anim.%s.gif", snn_str_i)
  ), anim_res)
  if (F) {
    # anime_res <- animate(panim+view_follow(), renderer = ffmpeg_renderer())
    # anim_save(file.path(dir_res,
    #                     sprintf('anim.%s.mp4', snn_str_i)), anim_res)
    anim_save(
      filename = file.path(dir_res, sprintf("anim.%s.mp4", snn_str_i)),
      animation = animate(panim + view_follow()),
      renderer = ffmpeg_renderer(format = "mp4", options = list(pix_fmt = "yuv420p", vcodec = "libx264"))
    )
  }
  gif_to_mp4 <- function(f_gif, to = NULL) {
    if (is.null(to)) {
      to <- paste0(f_gif, ".mp4")
    }
    if (file.exists(to)) {
      file_delete(to)
    }
    cmd <- paste0(
      "/usr/local/bin/ffmpeg ",
      "-i ", f_gif,
      " -movflags faststart ",
      " -pix_fmt yuv420p ",
      " -vf ",
      '\"scale=trunc(iw/2)*2:trunc(ih/2)*2\" ',
      to
    )
    cat(cmd)
    system(cmd)
  }
  try(gif_to_mp4(file.path(
    dir_res,
    sprintf("anim.%s.gif", snn_str_i)
  )))
}


animate_dimplot_spatial <- function(
    object,
    group.by = NULL,
    cells = NULL,
    cols = NULL) {

}

cat("Done\n")
