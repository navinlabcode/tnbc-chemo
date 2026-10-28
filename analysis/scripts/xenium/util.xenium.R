# Goal: helper functions for spatial
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yyan7@mdanderson.org)
# Date updated: Feb 3, 2026
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(igraph)
library(parallelDist)
library(ggrastr)
library(nanoparquet)
library(ggnewscale)

#------------------ ~~~ spatial visualization ~~~ --------------------

q_add_roi_box <- function(plot, df_rois, color = "black", linewidth = 1, labelsize = 1, labelpos = c("topleft", "center")) {
    labeladjust <- switch(match.arg(labelpos),
        topleft = c(0.1, 0.1),
        center = c(0.5, 0.5)
    )
    for (roi_index in seq_along(df_rois$ROI)) {
        roi_xmin <- df_rois$roi_xmin[roi_index]
        roi_xmax <- df_rois$roi_xmax[roi_index]
        roi_ymin <- df_rois$roi_ymin[roi_index]
        roi_ymax <- df_rois$roi_ymax[roi_index]
        roi_name <- as.character(df_rois$ROI[roi_index])
        plot <- plot + ggplot2::annotate("rect",
            xmin = roi_xmin, xmax = roi_xmax,
            ymin = roi_ymin, ymax = roi_ymax,
            color = color, fill = NA, linewidth = linewidth
        ) + ggplot2::annotate("text",
            x = roi_xmin + (roi_xmax - roi_xmin) * labeladjust[1],
            y = roi_ymin + (roi_ymax - roi_ymin) * labeladjust[2],
            label = roi_name,
            color = color, size = labelsize
        )
    }
    return(plot)
}


q_add_anno_text <- function(plot, text, color, roi_xmin, roi_xmax, roi_ymin, roi_ymax, size = 3) {
    plot + ggplot2::annotate("text",
        x = roi_xmin + (roi_xmax - roi_xmin) * c(0.9),
        y = roi_ymin + (roi_ymax - roi_ymin) * seq(from = 0.95, by = -0.05 / 3, length.out = length(text)),
        label = text,
        color = color, size = size
    )
}
q_add_anno_xenium_scale_bar <- function(plot, xum = NULL, roi_xmin = NULL, roi_xmax = NULL, roi_ymin = NULL, roi_ymax = NULL, color = "yellow", linewidth = 1) {
    if (is.null(xum)) {
        xyrange <- min(roi_xmax - roi_xmin, roi_ymax - roi_ymin)
        xum <- 10^floor(log10(xyrange))
    }
    plot + ggplot2::annotate("segment",
        x = roi_xmin + (roi_xmax - roi_xmin) * 0.05,
        xend = roi_xmin + (roi_xmax - roi_xmin) * 0.05 + xum,
        y = roi_ymin + (roi_ymax - roi_ymin) * 0.05,
        yend = roi_ymin + (roi_ymax - roi_ymin) * 0.05,
        color = color, linewidth = linewidth
    ) + ggplot2::annotate("text",
        x = roi_xmin + (roi_xmax - roi_xmin) * 0.05 + xum / 2,
        y = roi_ymin + (roi_ymax - roi_ymin) * 0.05 + xum * 0.1,
        label = paste0(xum, " um"),
        color = color, size = linewidth
    )
}

get_cellnames_ROI <- function(object, roi_xmin, roi_xmax, roi_ymin, roi_ymax) {
    df_coord <- GetTissueCoordinates(object)
    df_coord <- df_coord[df_coord$x >= roi_xmin & df_coord$x <= roi_xmax & df_coord$y >= roi_ymin & df_coord$y <= roi_ymax, ]
    return(df_coord$cell)
}

get_cell_segmentations_ROI <- function(object, roi_xmin, roi_xmax, roi_ymin, roi_ymax, fov = "fov") {
    cnames_roi <- get_cellnames_ROI(
        object = object,
        roi_xmin = roi_xmin, roi_xmax = roi_xmax,
        roi_ymin = roi_ymin, roi_ymax = roi_ymax
    )
    o <- FetchData(object[[fov]], vars = "segmentation", cells = cnames_roi)
    return(o)
}

write_xenium_rds <- function(object, dir_out, obj_type = c("pass", "ready")) {
    f_out <- switch(obj_type,
        pass = file.path(dir_out, "pass.xenium.seurat.rds"),
        ready = file.path(dir_out, "ready.xenium.seurat.rds")
    )
    f_meta_prefix <- switch(obj_type,
        pass = "pass.xenium.cellmeta.df",
        ready = "ready.xenium.cellmeta.df"
    )
    write_rds(object, f_out)
    cat("export object rds...")

    df <- cbind(
        GetTissueCoordinates(object) %>% column_to_rownames("cell"),
        object@meta.data
    )

    if ("umap" %in% names(object@reductions)) {
        df_red <- Embeddings(object, "umap")
        df <- cbind(df_red, df)
    }
    df$cellname <- Cells(object)

    write_rds(df, file.path(dir_out, paste0(f_meta_prefix, ".rds")))
    cat("export cell meta data rds...")
    write_csv(df, file.path(dir_out, paste0(f_meta_prefix, ".csv")))
    cat("export cell meta data csv...")
    nanoparquet::write_parquet(df, file.path(dir_out, paste0(f_meta_prefix, ".parquet")))
    cat("export cell meta data parquet...")

    cat("[done]\n")
}


ImageFeaturePlot.raster <- function(
    object, features, slot = c("data", "metadata"),
    size = 0.1, cex = 0.5,
    max.cutoff = "q95", min.cutoff = "q5",
    cols = c("lightgrey", "maroon1"),
    axes = FALSE, dark.background = FALSE,
    shuffle = 1026, order = FALSE) {
    if (length(features) > 1) {
        return(NULL)
    }
    slot <- match.arg(slot)
    df_coord <- GetTissueCoordinates(object)
    if (slot != "metadata") {
        df_coord$z <- deframe(FetchData(object, assay = DefaultAssay(object), slot = slot, vars = features))
    } else {
        df_coord$z <- object@meta.data[[features]]
    }


    if (as.numeric(shuffle) > 0) {
        set.seed(shuffle)
        df_coord <- df_coord[sample(1:nrow(df_coord), nrow(df_coord), replace = F), ]
    }
    if (order) {
        df_coord <- df_coord[order(df_coord$z), ]
    }

    if (!is.null(max.cutoff) & !is.null(min.cutoff)) {
        # exact value
        if (is.numeric(min.cutoff)) {
            z_min <- min.cutoff
        }
        if (is.numeric(max.cutoff)) {
            z_max <- max.cutoff
        }
        # text
        if (is.character(min.cutoff)) {
            min.cutoff <- min.cutoff %>%
                str_remove("q") %>%
                as.numeric()
            min.cutoff <- min.cutoff / 100
            z_min <- quantile(df_coord$z, probs = min.cutoff, na.rm = T)
        }
        if (is.character(max.cutoff)) {
            max.cutoff <- max.cutoff %>%
                str_remove("q") %>%
                as.numeric()
            max.cutoff <- max.cutoff / 100
            z_max <- quantile(df_coord$z, probs = max.cutoff, na.rm = T)
        }

        message(c(
            "Value range changes from ",
            min(df_coord$z), "-", max(df_coord$z), " to ",
            z_min, "-", z_max
        ))

        df_coord$z[df_coord$z > z_max] <- z_max
        df_coord$z[df_coord$z < z_min] <- z_min
    }

    p <- ggplot(data = df_coord, aes(x = x, y = y)) +
        geom_point_rast(aes(color = z), size = size, scale = cex) +
        scale_y_reverse() +
        coord_fixed() +
        labs(color = features)

    if (!axes) {
        p <- p + rremove("grid") + rremove("axis") + rremove("xylab") + rremove("xy.text") + rremove("ticks")
    }
    if (dark.background) {
        p <- p + theme(panel.background = element_rect(fill = "black", colour = "black"))
    }
    if (!all(is.null(cols))) {
        p <- p + scale_color_gradient(low = cols[1], high = cols[2])
    }
    return(p)
}

ImageDimPlot.raster <- function(
    object, group.by = NULL, cols = NULL,
    size = 0.1, cex = 0.5,
    roi_xmin = NULL, roi_xmax = NULL, roi_ymin = NULL, roi_ymax = NULL,
    axes = FALSE, dark.background = FALSE, shuffle = 1026,
    rev_yaxis_to_match_img = TRUE,
    add_roi_border = TRUE) {
    df_coord <- GetTissueCoordinates(object)
    df_coord$z <- object@meta.data[[group.by]]

    if (!"factor" %in% class(df_coord$z)) {
        df_coord$z <- as.factor(df_coord$z)
    }

    if (as.numeric(shuffle) > 0) {
        set.seed(shuffle)
        df_coord <- df_coord[sample(1:nrow(df_coord), nrow(df_coord), replace = F), ]
    } else {
        df_coord <- df_coord[rev(order(df_coord$z)), ]
    }

    if (!is.null(roi_xmin) & !is.null(roi_xmax) & !is.null(roi_ymin) & !is.null(roi_ymax)) {
        df_coord <- df_coord[df_coord$x >= roi_xmin & df_coord$x <= roi_xmax & df_coord$y >= roi_ymin & df_coord$y <= roi_ymax, ]
    }

    p <- ggplot(data = df_coord, aes(x = x, y = y)) +
        geom_point_rast(aes(color = z), size = size, scale = cex) +
        scale_color_manual(values = cols, na.value = cols["NA"]) +
        coord_fixed(expand = FALSE) +
        labs(color = group.by)
    if (add_roi_border) {
        p <- p + annotate("rect",
            ymin = roi_ymin,
            ymax = roi_ymax,
            xmin = roi_xmin,
            xmax = roi_xmax,
            color = ifelse(dark.background, "white", "black"), fill = NA
        )
    }
    if (!axes) {
        p <- p + rremove("grid") + rremove("axis") + rremove("xylab") + rremove("xy.text") + rremove("ticks")
    }
    if (dark.background) {
        p <- p + theme(panel.background = element_rect(fill = "black", colour = "black"))
    }
    if (rev_yaxis_to_match_img) {
        p <- p + scale_y_reverse()
    }
    return(p)
}


ImageDimPlotMol.raster <- function(
    object, molecules, mols.cols = NULL, nmols = Inf,
    fov = "fov",
    mols.size = 0.1, mols.cex = 0.5,
    color_mols_method = c("alpha_same", "alpha_data", "global_density", "local_density", "global_density2", "local_density2"),
    mols.alpha = NULL,
    roi_xmin = NULL, roi_xmax = NULL, roi_ymin = NULL, roi_ymax = NULL,
    add_cell_segmentation = FALSE, color_cell_segmentation_by = NULL,
    cell_segmentation_cols = NULL,
    cell_segmentation_width = 0.1,
    axes = FALSE, dark.background = FALSE, shuffle = 1026,
    rev_yaxis_to_match_img = TRUE,
    add_roi_border = TRUE, 
    do_mol_raster = TRUE) {
    
    func_geom_point <- geom_point_rast
    if (! do_mol_raster) {
        func_geom_point <- geom_point
    }

    df_coord <- FetchData(object[[fov]], vars = molecules, nmols = nmols)

    if (!"factor" %in% class(df_coord$molecule)) {
        df_coord$molecule <- as.factor(df_coord$molecule)
    }
    molecule_lvs <- levels(df_coord$molecule)
    if (is.null(mols.cols)) {
        mols.cols <- tidyplots::colors_discrete_candy[1:length(molecule_lvs)]
    }
    if (color_mols_method %in% c("global_density", "global_density2")) {
        df_coord <- split(df_coord, df_coord$molecule) %>%
            lapply(., function(x) {
                x$z <- get_density_val(x$x, x$y)
                return(x)
            }) %>%
            do.call(rbind, .)
    }

    if (!is.null(roi_xmin) & !is.null(roi_xmax) & !is.null(roi_ymin) & !is.null(roi_ymax)) {
        df_coord <- df_coord[df_coord$x >= roi_xmin & df_coord$x <= roi_xmax & df_coord$y >= roi_ymin & df_coord$y <= roi_ymax, ]
    }

    if (color_mols_method %in% c("local_density", "local_density2")) {
        # print(table(df_coord$molecule))
        df_coord <- split(df_coord, df_coord$molecule) %>%
            lapply(., function(x) {
                z <- try(get_density_val(x$x, x$y))
                if ("try-error" %in% class(z)) {
                    message("local density is fixed to 1")
                    z <- rep(1, nrow(x))
                }
                x$z <- z
                return(x)
            }) %>%
            do.call(rbind, .)
    }
    # print(head(df_coord))

    if (as.numeric(shuffle) > 0) {
        set.seed(shuffle)
        df_coord <- df_coord[sample(1:nrow(df_coord), nrow(df_coord), replace = F), ]
    } else {
        df_coord <- df_coord[rev(order(df_coord$molecule)), ]
    }

    df_seg <- NULL
    if (add_cell_segmentation) {
        df_seg <- get_cell_segmentations_ROI(
            object = object,
            roi_xmin = roi_xmin, roi_xmax = roi_xmax,
            roi_ymin = roi_ymin, roi_ymax = roi_ymax, fov = fov
        )
    }
    if (!is.null(color_cell_segmentation_by)) {
        df_cell_groupby <- FetchData(object, vars = color_cell_segmentation_by)
        df_cell_groupby$z <- df_cell_groupby[[color_cell_segmentation_by]]
        df_cell_groupby[[color_cell_segmentation_by]] <- NULL
        df_cell_groupby$cell <- rownames(df_cell_groupby)
        df_seg <- left_join(df_seg, df_cell_groupby, by = "cell")
    }
    p <- ggplot()
    if (!is.null(df_seg)) {
        if (!is.null(color_cell_segmentation_by)) {
            p <- p + geom_polygon(
                data = df_seg,
                aes(x = x, y = y, group = cell, fill = z),
                linewidth = cell_segmentation_width,
                color = ifelse(dark.background, "black", "white")
            ) +
                scale_fill_manual(values = cell_segmentation_cols) +
                ggnewscale::new_scale_color()
        } else {
            p <- p + geom_polygon(
                data = df_seg,
                aes(x = x, y = y, group = cell),
                linewidth = cell_segmentation_width, fill = NA,
                color = ifelse(dark.background, "white", "black")
            )
        }
    }
    if (color_mols_method == "alpha_same") {
        if (is.null(mols.alpha)) {
            mols.alpha <- 0.5
        }
        p <- p +
            func_geom_point(
                data = df_coord,
                aes(x = x, y = y, color = molecule),
                size = mols.size, pch = 16,
                scale = mols.cex, alpha = mols.alpha
            ) +
            scale_color_manual(values = mols.cols)
    }
    if (color_mols_method == "alpha_data") {
        # todo
    }
    if (color_mols_method %in% c("global_density", "local_density")) {
        for (mm in molecule_lvs) {
            p <- p +
                func_geom_point(
                    data = df_coord[df_coord$molecule == mm, ],
                    aes(x = x, y = y, color = z),
                    size = mols.size, pch = 16,
                    scale = mols.cex
                ) +
                scale_color_gradient(
                    low = ifelse(dark.background, "grey50", "grey90"),
                    high = mols.cols[mm], name = mm
                ) +
                ggnewscale::new_scale_color()
        }
    }
    # print(head(df_coord))
    if (color_mols_method %in% c("global_density2", "local_density2")) {
        for (mm in molecule_lvs) {
            p <- p +
                func_geom_point(
                    data = df_coord[df_coord$molecule == mm, ],
                    aes(x = x, y = y, alpha = z),
                    size = mols.size, pch = 16,
                    scale = mols.cex, color = mols.cols[mm]
                ) +
                ggnewscale::new_scale_color()
        }
    }
    if (add_roi_border) {
        p <- p + annotate("rect",
            ymin = roi_ymin,
            ymax = roi_ymax,
            xmin = roi_xmin,
            xmax = roi_xmax,
            color = ifelse(dark.background, "white", "black"), fill = NA
        )
    }
    if (!axes) {
        p <- p + rremove("grid") + rremove("axis") + rremove("xylab") + rremove("xy.text") + rremove("ticks")
    }
    if (dark.background) {
        p <- p + theme(panel.background = element_rect(fill = "black", colour = "black"))
    }
    if (rev_yaxis_to_match_img) {
        p <- p +
            scale_y_reverse() +
            coord_fixed(
                expand = FALSE,
                xlim = c(roi_xmin, roi_xmax),
                ylim = c(roi_ymax, roi_ymin)
            )
    } else {
        p <- p + coord_fixed(
            expand = FALSE,
            xlim = c(roi_xmin, roi_xmax),
            ylim = c(roi_ymin, roi_ymax)
        )
    }
    return(p)
}


ImageDimPlotSegmentation.raster <- function(
    object, group.by = NULL, cols = NULL, fov = "fov",
    size = 0.1, cex = 0.5, alpha = 1,
    cell_segmentation_width = 0.01,
    roi_xmin = NULL, roi_xmax = NULL, roi_ymin = NULL, roi_ymax = NULL,
    axes = FALSE, dark.background = FALSE, shuffle = 1026,
    rev_yaxis_to_match_img = TRUE,
    add_roi_border = TRUE, do_raster = TRUE) {
    df_sp <- GetTissueCoordinates(object)
    if (is.null(roi_xmin)) {
        roi_xmin <- min(df_sp$x)
    }
    if (is.null(roi_xmax)) {
        roi_xmax <- max(df_sp$x)
    }
    if (is.null(roi_ymin)) {
        roi_ymin <- min(df_sp$y)
    }
    if (is.null(roi_ymax)) {
        roi_ymax <- max(df_sp$y)
    }

    df_coord <- get_cell_segmentations_ROI(
        object = object,
        roi_xmin = roi_xmin, roi_xmax = roi_xmax,
        roi_ymin = roi_ymin, roi_ymax = roi_ymax, fov = fov
    ) %>% as_tibble()
    df_coord$cell <- as.character(df_coord$cell)
    cnames <- unique(df_coord$cell)

    df_coord_z <- FetchData(object, vars = group.by, cells = cnames)
    df_coord_z$cell <- rownames(df_coord_z)

    df_coord <- left_join(df_coord, df_coord_z, by = "cell")
    df_coord$z <- df_coord[[group.by]]

    if (!"factor" %in% class(df_coord$z)) {
        df_coord$z <- as.factor(df_coord$z)
    }
    if (as.numeric(shuffle) > 0) {
        set.seed(shuffle)
        df_coord$cell <- factor(df_coord$cell, levels = sample(cnames, size = length(cnames), replace = F))
    } else {
        df_coord <- df_coord %>% dplyr::arrange(z, cell)
    }

    if (do_raster) {
        # options(ggrastr.default.dpi=1000)
        p <- ggplot(data = df_coord, aes(x = x, y = y)) +
            ggrastr::rasterise(geom_polygon(
                aes(group = cell, fill = z),
                linewidth = cell_segmentation_width, alpha = alpha,
                color = ifelse(dark.background, "black", "white")
            )) +
            scale_fill_manual(values = cols, na.value = cols["NA"]) +
            labs(fill = group.by)
    } else {
        p <- ggplot(data = df_coord, aes(x = x, y = y)) +
            geom_polygon(
                aes(group = cell, fill = z),
                linewidth = cell_segmentation_width, alpha = alpha,
                color = ifelse(dark.background, "black", "white")
            ) +
            scale_fill_manual(values = cols, na.value = cols["NA"]) +
            labs(fill = group.by)
    }


    if (add_roi_border) {
        p <- p + annotate("rect",
            ymin = roi_ymin,
            ymax = roi_ymax,
            xmin = roi_xmin,
            xmax = roi_xmax,
            color = ifelse(dark.background, "white", "black"), fill = NA
        )
    }
    if (!axes) {
        p <- p + rremove("grid") + rremove("axis") + rremove("xylab") + rremove("xy.text") + rremove("ticks")
    }
    if (dark.background) {
        p <- p + theme(panel.background = element_rect(fill = "black", colour = "black"))
    }
    if (rev_yaxis_to_match_img) {
        p <- p + scale_y_reverse() +
            coord_fixed(expand = FALSE, xlim = c(roi_xmin, roi_xmax), ylim = c(roi_ymax, roi_ymin))
    } else {
        p <- p +
            coord_fixed(expand = FALSE, xlim = c(roi_xmin, roi_xmax), ylim = c(roi_ymin, roi_ymax))
    }
    return(p)
}


ImageDimPlotSegmentation.raster.2 <- function(
    object, group.by = NULL, cols = NULL, fov = "fov",
    # size = 0.1, cex = 0.5,
    alpha = 1,
    cell_segmentation_by = NULL,
    cell_segmentation_cols = NULL,
    cell_segmentation_width = 0.1,
    roi_xmin = NULL, roi_xmax = NULL, roi_ymin = NULL, roi_ymax = NULL,
    axes = FALSE, dark.background = FALSE, shuffle = 1026,
    rev_yaxis_to_match_img = TRUE,
    add_roi_border = TRUE, do_raster = TRUE) {
    df_sp <- GetTissueCoordinates(object)
    if (is.null(roi_xmin)) {
        roi_xmin <- min(df_sp$x)
    }
    if (is.null(roi_xmax)) {
        roi_xmax <- max(df_sp$x)
    }
    if (is.null(roi_ymin)) {
        roi_ymin <- min(df_sp$y)
    }
    if (is.null(roi_ymax)) {
        roi_ymax <- max(df_sp$y)
    }

    df_coord <- get_cell_segmentations_ROI(
        object = object,
        roi_xmin = roi_xmin, roi_xmax = roi_xmax,
        roi_ymin = roi_ymin, roi_ymax = roi_ymax, fov = fov
    ) %>% as_tibble()
    df_coord$cell <- as.character(df_coord$cell)
    cnames <- unique(df_coord$cell)

    df_coord_z <- FetchData(object, vars = group.by, cells = cnames)
    df_coord_z$cell <- rownames(df_coord_z)
    df_coord_e <- FetchData(object, vars = cell_segmentation_by, cells = cnames)
    df_coord_e$cell <- rownames(df_coord_e)

    df_coord <- left_join(df_coord, df_coord_z, by = "cell")
    df_coord <- left_join(df_coord, df_coord_e, by = "cell")
    df_coord$z <- df_coord[[group.by]]
    df_coord$e <- df_coord[[cell_segmentation_by]]

    if (!"factor" %in% class(df_coord$z)) {
        df_coord$z <- as.factor(df_coord$z)
    }
    if (!"factor" %in% class(df_coord$e)) {
        df_coord$e <- as.factor(df_coord$e)
    }
    # print(table(df_coord$e))
    if (as.numeric(shuffle) > 0) {
        set.seed(shuffle)
        df_coord$cell <- factor(df_coord$cell, levels = sample(cnames, size = length(cnames), replace = F))
    } else {
        df_coord <- df_coord %>% dplyr::arrange(z, cell)
    }
    if (do_raster) {
        p <- ggplot(data = df_coord, aes(x = x, y = y)) +
            ggrastr::rasterise(
                geom_polygon(
                    aes(group = cell, fill = z),
                    color = NA,
                    linewidth = cell_segmentation_width, alpha = alpha
                ),
                dev = "cairo"
            ) +
            scale_fill_manual(values = cols, na.value = cols["NA"]) +
            ggrastr::rasterise(
                geom_polygon(
                    aes(group = cell, color = e),
                    fill = NA,
                    linewidth = cell_segmentation_width, alpha = alpha
                ),
                dev = "cairo"
            ) +
            scale_color_manual(values = cell_segmentation_cols, na.value = cell_segmentation_cols["NA"]) +
            labs(fill = group.by, color = cell_segmentation_by)
    } else {
        p <- ggplot(data = df_coord, aes(x = x, y = y)) +
            geom_polygon(
                aes(group = cell, fill = z),
                color = NA,
                linewidth = cell_segmentation_width, alpha = alpha
            ) +
            scale_fill_manual(values = cols, na.value = cols["NA"]) +
            geom_polygon(
                aes(group = cell, color = e),
                fill = NA,
                linewidth = cell_segmentation_width, alpha = alpha
            ) +
            scale_color_manual(values = cell_segmentation_cols, na.value = cell_segmentation_cols["NA"]) +
            labs(fill = group.by, color = cell_segmentation_by)
    }

    if (add_roi_border) {
        p <- p + annotate("rect",
            ymin = roi_ymin,
            ymax = roi_ymax,
            xmin = roi_xmin,
            xmax = roi_xmax,
            color = ifelse(dark.background, "white", "black"), fill = NA
        )
    }
    if (!axes) {
        p <- p + rremove("grid") + rremove("axis") + rremove("xylab") + rremove("xy.text") + rremove("ticks")
    }
    if (dark.background) {
        p <- p + theme(panel.background = element_rect(fill = "black", colour = "black"))
    }
    if (rev_yaxis_to_match_img) {
        p <- p + scale_y_reverse() +
            coord_fixed(expand = FALSE, xlim = c(roi_xmin, roi_xmax), ylim = c(roi_ymax, roi_ymin))
    } else {
        p <- p +
            coord_fixed(expand = FALSE, xlim = c(roi_xmin, roi_xmax), ylim = c(roi_ymin, roi_ymax))
    }
    return(p)
}

fetch_spatial_coord <- function(object, cells, fov) {
    if (missing(fov)) {
        fov <- "fov"
    }
    if (missing(cells)) {
        cells <- Cells(object)
    }
    fov <- object[[fov]]
    df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
    df_coord <- df_coord[cells, ]
    return(df_coord)
}
get_spatial_xy_ratio <- function(object, cells, fov) {
    if (missing(fov)) {
        fov <- "fov"
    }
    if (missing(cells)) {
        cells <- Cells(object)
    }
    fov <- object[[fov]]
    df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
    df_coord <- df_coord[cells, ]
    coord_xy_ratio <- diff(range(df_coord$x)) / diff(range(df_coord$y))
    return(coord_xy_ratio)
}


to_img_xenium_pdf_inch <- function(x) {
    x / 1e2
}

calc_img_spatial_pdf_size <- function(spatial_xy_ratio, fixed_inch = 7.5, flip_xy = FALSE) {
    if (flip_xy) {
        spatial_xy_ratio <- 1 / spatial_xy_ratio
    }
    if (spatial_xy_ratio < 1) {
        pdf_img_height <- fixed_inch
        pdf_img_width <- pdf_img_height * spatial_xy_ratio
        pdf_img_width <- pmax(pdf_img_width, 3)
    } else {
        pdf_img_width <- fixed_inch
        pdf_img_height <- pdf_img_width / spatial_xy_ratio
        pdf_img_height <- pmax(pdf_img_height, 3)
    }
    return(c("height" = pdf_img_height, "width" = pdf_img_width))
}

quick_spatial_distance <- function(object, cells, fov, ...) {
    if (missing(cells)) {
        cells <- Cells(object)
    }
    if (missing(fov)) {
        fov <- "fov"
    }
    df_coord <- fetch_spatial_coord(object, cells = cells)
    parallelDist::parallelDist(as.matrix(df_coord), ...)
}

make_seurat_flip_xy_FALSE_same_as_in_situ <- scale_x_reverse()
make_raw_spatial_same_as_in_situ <- scale_y_reverse()

ImgDimPlotAnim <- function(object, group.by, to_gif_path, cols = NULL, do_downsample = F) {
    Idents(object) <- group.by
    fov <- object[["fov"]]
    head(fov$centroids@coords)
    head(fov$centroids@cells)
    df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
    df_coord <- df_coord[Cells(object), ]
    coord_xy_ratio <- diff(range(df_coord$x)) / diff(range(df_coord$y))
    df_coord$x <- rescale(df_coord$x, to = c(-1 * coord_xy_ratio, coord_xy_ratio))
    df_coord$y <- -1 * rescale(df_coord$y, to = c(-1, 1))
    df_spatial <- data.frame(df_coord,
        cluster = Idents(object),
        type = "spatial", stringsAsFactors = F
    )

    df_umap <- as.data.frame(Embeddings(object, "umap"))
    colnames(df_umap) <- c("x", "y")
    df_umap$x <- rescale(df_umap$x, to = c(-1, 1))
    df_umap$y <- rescale(df_umap$y, to = c(-1, 1))
    df_umap <- data.frame(df_umap,
        cluster = Idents(object),
        type = "umap", stringsAsFactors = F
    )

    if (do_downsample) {
        tmp <- sample(1:nrow(df_spatial), size = round(nrow(df_spatial) / 3), replace = F)
        df_trans <- rbind(df_spatial[tmp, ], df_umap[tmp, ])
    } else {
        df_trans <- rbind(df_spatial, df_umap)
    }

    if (is.null(cols)) {
        cols <- structure(
            Seurat::DiscretePalette(n = length(levels(Idents(object))), palette = "parade"),
            names = levels(Idents(object))
        )
    }
    panim <- ggplot(df_trans, aes(x = x, y = y)) +
        geom_scattermore(
            aes(color = seurat_clusters),
            pointsize = 1.2, pixels = c(1024, 1024)
        ) +
        scale_color_manual(values = cols) +
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
    anim_save(to_gif_path, anim_res)
    return(invisible(anim_res))
}

ImgFeaturePlotAnim <- function(object, features, to_gif_path, cols = NULL,
                               slot = "data",
                               do_downsample = F) {
    fov <- object[["fov"]]
    head(fov$centroids@coords)
    head(fov$centroids@cells)

    z <- deframe(FetchData(object, slot = slot, vars = features))

    df_coord <- data.frame(fov$centroids@coords, row.names = fov$centroids@cells)
    df_coord <- df_coord[Cells(object), ]
    coord_xy_ratio <- diff(range(df_coord$x)) / diff(range(df_coord$y))
    df_coord$x <- rescale(df_coord$x, to = c(-1 * coord_xy_ratio, coord_xy_ratio))
    df_coord$y <- -1 * rescale(df_coord$y, to = c(-1, 1))
    df_spatial <- data.frame(
        df_coord,
        val = z,
        type = "spatial", stringsAsFactors = F
    )

    df_umap <- as.data.frame(Embeddings(object, "umap"))
    colnames(df_umap) <- c("x", "y")
    df_umap$x <- rescale(df_umap$x, to = c(-1, 1))
    df_umap$y <- rescale(df_umap$y, to = c(-1, 1))
    df_umap <- data.frame(
        df_umap,
        val = z,
        type = "umap", stringsAsFactors = F
    )

    if (do_downsample) {
        tmp <- sample(1:nrow(df_spatial), size = round(nrow(df_spatial) / 3), replace = F)
        df_trans <- rbind(df_spatial[tmp, ], df_umap[tmp, ])
    } else {
        df_trans <- rbind(df_spatial, df_umap)
    }

    panim <- ggplot(df_trans, aes(x = x, y = y)) +
        geom_scattermore(
            aes(x = x, y = y, color = z),
            pointsize = 1.2, pixels = c(1024, 1024)
        ) +
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
    if (!all(is.null(cols))) {
        panim <- panim + scale_color_gradient(low = cols[1], high = cols[2])
    }
    panim <- panim + rremove("legend")

    panim <- panim +
        transition_states(type,
            transition_length = 5,
            state_length = 1
        ) +
        labs(title = "{closest_state}") +
        theme(plot.title = element_text(size = 28)) +
        enter_fade()

    anim_res <- panim + view_follow()
    anim_save(to_gif_path, anim_res)
    return(invisible(anim_res))
}

#------------------ ~~~ Xenium Explorer ~~~ --------------------
df_to_exnium_explorer <- function(df, cell_str, cat_str) {
    df_new <- df[, c(cell_str, cat_str)]
    colnames(df_new) <- c("cell_id", "group")
    return(df_new)
}

xenium_explorer_selection_to_dataframe <- function(x, check_file = FALSE) {
    df <- read_csv(x, comment = "#", skip = 2, show_col_types = FALSE)
    if (check_file) {
        names_rois_requested <- read_lines(x, n_max = 1) %>%
            str_remove("#Selection names:,") %>%
            str_split_1(",") %>%
            str_trim() %>%
            str_remove_all("\"")
        if (!identical(sort(unique(df$Selection)), sort(names_rois_requested))) {
            warning(
                "The requested ROIs are not identical to the ROIs in the file. Please check the file."
            )
        }
    }
    df_wide <- df %>%
        group_by(Selection) %>%
        summarise(roi_xmin = min(X), roi_xmax = max(X), roi_ymin = min(Y), roi_ymax = max(Y)) %>%
        ungroup()
    df_wide$ROI <- str_replace_all(df_wide$Selection, "Selection", "ROI") %>%
        str_squish() %>%
        str_replace_all(" ", "_")
    names_rois_requested <- df_wide$ROI
    df_wide$ROI <- factor(df_wide$ROI, levels = gtools::mixedsort(names_rois_requested))
    df_wide <- df_wide[order(df_wide$ROI), ]
    return(df_wide)
}

#------ Helper funcs for Xenium color bar ------
q_add_anno_text <- function(plot, text, color, roi_xmin, roi_xmax, roi_ymin, roi_ymax, size = 3) {
    plot + ggplot2::annotate("text",
        x = roi_xmin + (roi_xmax - roi_xmin) * c(0.9),
        y = roi_ymin + (roi_ymax - roi_ymin) * seq(from = 0.95, by = -0.05 / 3, length.out = length(text)),
        label = text,
        color = color, size = size
    )
}
q_add_anno_xenium_scale_bar <- function(plot, xum = NULL, roi_xmin = NULL, roi_xmax = NULL, roi_ymin = NULL, roi_ymax = NULL, color = "yellow", linewidth = 1) {
    if (is.null(xum)) {
        xyrange <- min(roi_xmax - roi_xmin, roi_ymax - roi_ymin)
        xum <- 10^floor(log10(xyrange))
    }
    plot + ggplot2::annotate("segment",
        x = roi_xmin + (roi_xmax - roi_xmin) * 0.05,
        xend = roi_xmin + (roi_xmax - roi_xmin) * 0.05 + xum,
        y = roi_ymin + (roi_ymax - roi_ymin) * 0.05,
        yend = roi_ymin + (roi_ymax - roi_ymin) * 0.05,
        color = color, linewidth = linewidth
    ) + ggplot2::annotate("text",
        x = roi_xmin + (roi_xmax - roi_xmin) * 0.05 + xum / 2,
        y = roi_ymin + (roi_ymax - roi_ymin) * 0.05 + xum * 0.1,
        label = paste0(xum, " um"),
        color = color, size = linewidth
    )
}
