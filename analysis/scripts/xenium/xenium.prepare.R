suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: given a xenium object, perform the lazy preprocessing
# normalization, scaledata, runpca, runumap, etc...
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-08-17
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
    library(ruok)
    library(Seurat)
    # library(arrow)
    library(fs)
    library(Signac)
    library(clustree)
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
    my_scatter_themevoid <- theme_pubr(base_size = 6, legend = "right") %+replace% theme(
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = rel(1)),
        axis.line = element_blank()
    )
})
options <- commandArgs(trailingOnly = TRUE)

#------------------ ~~~ input ~~~ --------------------

if (length(options) == 1) {
    f_in <- options[[1]]
    dir_res <- dirname(f_in) # file.path(, 'lognorm')
    assay_use <- "Xenium"
}
if (length(options) == 2) {
    f_in <- options[[1]]
    dir_res <- options[[2]]
    assay_use <- "Xenium"
}
if (length(options) == 3) {
    f_in <- options[[1]]
    dir_res <- options[[2]]
    assay_use <- options[[3]]
}


n_linear_dr_components <- 50
n_linear_dr_components_for_umap <- 50

#------ find input and setup output ------

f_out <- file.path(dir_res, "ready.seurat.rds")

dir_create(dir_res)
message(dir_res)
message(f_out)

if (file_exists(f_out)) {
# if (F) {
    message("Load the existing object...")
    xmo <- read_rds(f_out)
    DefaultAssay(xmo) <- assay_use
    print(xmo)
} else {
    message("Processing begins...")
    #------------------ ~~~ Preprocessing ~~~ --------------------
    xmo <- read_rds(f_in)
    print(xmo)
    DefaultAssay(xmo) <- assay_use

    if (assay_use != "integrated") {
        xmo <- NormalizeData(xmo)
    }

    xmo <- FindVariableFeatures(xmo)

    is_non_empty_features <- rowSums(GetAssayData(xmo, layer = "counts", assay = "Xenium")) != 0
    non_empty_features <- VariableFeatures(xmo) <- rownames(GetAssayData(xmo, layer = "counts", assay = "Xenium"))[is_non_empty_features]
    print(xmo)
    print(DefaultAssay(xmo))

    xmo <- ScaleData(
        xmo,
        vars.to.regress = intersect(sprintf("nCount_%s", "Xenium"), colnames(xmo@meta.data))
    )
    xmo <- RunPCA(xmo, npcs = n_linear_dr_components)
    write_rds(xmo, f_out)

    #------ umap ------
    i_linear_dr_components_for_umap <- 1:n_linear_dr_components_for_umap
    xmo <- RunUMAP(xmo,
        reduction = "pca",
        dims = i_linear_dr_components_for_umap
    )
    write_seurat(xmo, dir_res, obj_type = "ready")
    #------ cell clustering ------
    xmo <- FindNeighbors(
        xmo,
        k.param = 20,
        dims = i_linear_dr_components_for_umap,
        reduction = "pca"
    )
    write_seurat(xmo, dir_res, obj_type = "ready")

    idx <- grepl(pattern = "_snn_res", x = colnames(xmo[[]]))
    for (i in colnames(xmo[[]])[idx]) {
        xmo[[i]] <- NULL
    }
    snn_res_max <- 0.8
    if (ncol(xmo) < 100) {
        snn_res_max <- .2
    }

    for (snn_res_i in seq(from = 0.2, to = snn_res_max, by = .2)) {
        cat(snn_res_i, "... ")
        try(xmo <- FindClusters(
            xmo,
            algorithm = 3,
            resolution = snn_res_i, verbose = F
        ))
        snn_str_i <- sprintf("%s_snn_res.%s", DefaultAssay(xmo), snn_res_i)
        if (snn_str_i %in% colnames(xmo[[]])) {
            xmo[[snn_str_i]] <- Idents(xmo) ## The cluster order is human-readable
        }
    }
    cat("\n")
    write_seurat(xmo, dir_res, obj_type = "ready")

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
if ("FOV" %in% class(try(xmo[["fov"]]))) {
    pdf_img_width <- calc_img_spatial_pdf_size(get_spatial_xy_ratio(xmo))["width"]
    pdf_img_height <- calc_img_spatial_pdf_size(get_spatial_xy_ratio(xmo))["height"]
}

#------ viz technical features ------

for (feature_i in c(
    "nCount_Xenium", "nFeature_Xenium",
    "nCount_SCT", "nFeature_SCT",
    "nCount_integrated", "nFeature_integrated",
    "nCount_RNA", "nFeature_RNA",
    "EPCAM"
)) {
    if (!(feature_i %in% colnames(xmo@meta.data) | feature_i %in% rownames(xmo))) {
        cat("requested feature ", feature_i, "is not available...\n")
        next()
    }
    message(feature_i)

    dir_snippet_viz <- dir_res
    print(xmo)
    pal_v <- c("lightgrey", "blue") ## either c('low_color', 'high_color') or just NULL (to be change outside)
    viz_which_val <- feature_i

    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_continuous.R")
}

#------ viz idents ------

snn_res_i <- 0.2
snn_str_i <- sprintf("%s_snn_res.%s", DefaultAssay(xmo), snn_res_i)

snn_str_opts <- colnames(xmo@meta.data)[grepl(sprintf("%s_snn_res", DefaultAssay(xmo)), colnames(xmo@meta.data))]
viz_z_opts <- c(snn_str_opts, c("sample"))
viz_z_opts <- intersect(viz_z_opts, colnames(xmo@meta.data))
for (snn_str_i in viz_z_opts) {
    xmo$seurat_clusters <- xmo@meta.data[, snn_str_i]

    if (!"factor" %in% class(xmo$seurat_clusters)) {
        xmo$seurat_clusters <- factor(
            xmo$seurat_clusters,
            levels = gtools::mixedsort(unique(as.character(xmo$seurat_clusters)))
        )
    }

    pal_snn <- structure(
        Seurat::DiscretePalette(n = length(levels(xmo$seurat_clusters)), palette = "parade"),
        names = levels(xmo$seurat_clusters)
    )

    dir_snippet_viz <- dir_res
    viz_what <- snn_str_i
    pal_z <- pal_snn
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_categorical.R")
}

#------ animation of UMAP to spacial coordinates ------
#
# if (! 'FOV' %in% class(try(xmo[['fov']]))) {
#   viz_z_opts <- c()
# }
# library(scattermore)
# library(gganimate)
# for (snn_str_i in viz_z_opts) {
#   xmo$seurat_clusters <- xmo@meta.data[, snn_str_i]
#
#   pal_snn <- structure(
#     Seurat::DiscretePalette(n=length(levels(xmo$seurat_clusters)), palette = 'parade'),
#     names=levels(xmo$seurat_clusters))
#
#   fov <- xmo[['fov']]
#   head(fov$centroids@coords)
#   head(fov$centroids@cells)
#   df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
#   df_coord <- df_coord[Cells(xmo), ]; coord_xy_ratio <- diff(range(df_coord$x)) / diff(range(df_coord$y))
#   df_coord$x <- rescale(df_coord$x, to=c(-1*coord_xy_ratio, coord_xy_ratio))
#   df_coord$y <- -1 * rescale(df_coord$y, to=c(-1, 1))
#   df_spatial <- data.frame(df_coord, seurat_clusters=xmo$seurat_clusters, type='spatial', stringsAsFactors = F)
#
#   df_umap <- as.data.frame(Embeddings(xmo, 'umap'))
#   colnames(df_umap) <- c('x', 'y')
#   df_umap$x <- rescale(df_umap$x, to=c(-1, 1))
#   df_umap$y <- rescale(df_umap$y, to=c(-1, 1))
#   df_umap <- data.frame(df_umap, seurat_clusters=xmo$seurat_clusters, type='embedding', stringsAsFactors = F)
#
#   tmp <- sample(1:nrow(df_spatial), size = round(nrow(df_spatial)/3), replace = F)
#   df_trans <- rbind(df_spatial[tmp, ], df_umap[tmp, ])
#   # df_trans <- rbind(df_spatial, df_umap)
#
#   panim <- ggplot(df_trans, aes(x=x, y=y)) +
#     # geom_scattermore(
#     #   aes(x=x,y=y, color=seurat_clusters), pointsize=.5) +
#     geom_scattermore(
#       aes(x=x,y=y, color=seurat_clusters), pointsize=1.2, pixels = c(1024,1024)) +
#     scale_color_manual(values = pal_snn) +
#     Seurat::DarkTheme() +
#     rremove('x.title') + rremove('y.title') +
#     rremove('x.text') + rremove('y.text') + rremove('x.axis') + rremove('y.axis') +
#     coord_equal() +
#     theme(legend.position = "none")
#   # panim
#
#   panim <- panim +
#     transition_states(type,
#                       transition_length = 5,
#                       state_length = 1) +
#     labs(title = '{closest_state}') +
#     theme(plot.title = element_text(size = 28)) +
#     enter_fade()
#
#   anim_res <- panim + view_follow()
#   anim_save(file.path(dir_res,
#                       sprintf('anim.%s.gif', snn_str_i)), anim_res)
#   if(F){
#     # anime_res <- animate(panim+view_follow(), renderer = ffmpeg_renderer())
#     # anim_save(file.path(dir_res,
#     #                     sprintf('anim.%s.mp4', snn_str_i)), anim_res)
#     anim_save(filename = file.path(dir_res, sprintf('anim.%s.mp4', snn_str_i)),
#               animation =  animate(panim+view_follow()),
#               renderer = ffmpeg_renderer(format = "mp4", options = list(pix_fmt = "yuv420p", vcodec="libx264")))
#
#   }
#   gif_to_mp4 <- function(f_gif, to=NULL) {
#     if (is.null(to))  { to <- paste0(f_gif, '.mp4') }
#     if (file.exists(to)) { file_delete(to) }
#     cmd = paste0('/usr/local/bin/ffmpeg ',
#                  '-i ', f_gif,
#                  ' -movflags faststart ',
#                  ' -pix_fmt yuv420p ',
#                  ' -vf ' ,
#                  '\"scale=trunc(iw/2)*2:trunc(ih/2)*2\" ',
#                  to)
#     cat(cmd)
#     system(cmd)
#
#   }
#   try(gif_to_mp4(file.path(dir_res,
#                            sprintf('anim.%s.gif', snn_str_i))))
# }

cat("DONE xenium.prepare.R\n")
