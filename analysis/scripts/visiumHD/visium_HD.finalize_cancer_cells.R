suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal:Use RCTD to identify cell types in Visium HD data
# Create a data frame of cell types depending on the reference library
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Chenling Tang; Yun Yan (yun.yan@uth.tmc.edu)
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
    library(spacexr)
    library(Rfast)
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.visiumHD.R")
})
cmdargs <- commandArgs(trailingOnly = TRUE)

if (length(cmdargs) > 0) {
    f_in <- cmdargs[1]
    sample_name <- cmdargs[2]
} else {
    sample_name <- "ART122"
    sample_name <- "ART266"
    f_in <- glue("/volumes/USR1/yyan/project/tnbc_visium_hd/data/{sample_name}/binsize_32/ready.seurat.rds")
}

#------------------ ~~~ REadin ~~~ --------------------
cli_h1("REadin")

sr3 <- read_rds(f_in)
print(sr3)

dir_res <- file.path(dirname(f_in), sprintf("finalize_cancer_cells"))
fs::dir_create(dir_res)

#------------------ ~~~ Collect RCTD based on HBCA and TNBC ref ~~~ --------------------
cli_h1("Collect RCTD based on HBCA and TNBC ref")
rctd_res_df_hbca <- file.path(dirname(f_in), "celltype_infer_from_HBCA", "RCTD.seurat_dataframe.rds") %>% read_rds()
rctd_res_df_tnbc <- file.path(dirname(f_in), "celltype_infer_from_TNBC", "RCTD.seurat_dataframe.rds") %>% read_rds()

stopifnot(all(rownames(rctd_res_df_hbca) == rownames(rctd_res_df_tnbc)))
stopifnot(all(rownames(rctd_res_df_hbca) == Cells(sr3)))

#------------------ ~~~ Raw CopyKat results ~~~ --------------------
## should have more cells than sr3
cli_h1("CopyKat results")
copykat_res_df <- file.path(dirname(f_in), "copykat_mix", "copykat_pred_report.seurat3_meta.rds") %>% read_rds()
str(rownames(copykat_res_df))
str(Cells(sr3))
print(colnames(copykat_res_df))
shared_ncells <- intersect(Cells(sr3), rownames(copykat_res_df))
if (!identical(Cells(sr3), rownames(copykat_res_df))) {
    cli::cli_alert_warning("[expected] Mismatch in cell names between CopyKat and Seurat object.")
    tmp <- intersect(Cells(sr3), rownames(copykat_res_df))
    if (length(tmp) / length(Cells(sr3)) < 0.1) {
        stop("Less than 10% of cells in Seurat object are found in CopyKat results. Something must be wrong [error]")
    }

    idx <- match(Cells(sr3), rownames(copykat_res_df))
    copykat_res_df <- copykat_res_df[idx, , drop = FALSE]
    rownames(copykat_res_df) <- Cells(sr3)

    colnames(copykat_res_df)
    copykat_res_df$copykat_pred_default <- replace_na(copykat_res_df$copykat_pred_default, "Unknown")
    copykat_res_df$copykat_pred_tirosh <- replace_na(copykat_res_df$copykat_pred_tirosh, "Unknown")
    copykat_res_df$copykat_pred_leiden <- replace_na(copykat_res_df$copykat_pred_leiden, "Unknown")
    copykat_res_df$icna_leiden <- replace_na(as.character(copykat_res_df$icna_leiden), "Unknown") %>% as.factor()
}
stopifnot(all(Cells(sr3) == rownames(copykat_res_df)))

write_rds(copykat_res_df, file.path(dir_res, "copykat_pred_report.seurat3_meta.rds"))
write.csv(copykat_res_df, file.path(dir_res, "copykat_pred_report.seurat3_meta.csv"))

#------------------ ~~~ Viz ~~~ --------------------
cli_h1("Viz")
library(ruok)
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
## table the cell types of HBCA, TNBC and CopyKat

#------------------ ~~~ Comparing results ~~~ --------------------
cli_h2("Comparing results")
#------ TNBC vs HBCA ------
p1 <- ruok::qtable_scatter(rctd_res_df_tnbc$first_type, rctd_res_df_hbca$first_type) +
    labs(x = "TNBC", y = "HBCA") +
    scale_color_manual(values = pal_celltypes)
ggsave(file.path(dir_res, "scatter.celltypes_TNBC_vs_HBCA.pdf"), p1, width = 6, height = 5, useDingbats = FALSE)

p2 <- ruok::qbarplot_table_catx(rctd_res_df_tnbc$first_type, rctd_res_df_hbca$first_type, name_h = "TNBC cell types", name_v = "HBCA") +
    scale_fill_manual(values = pal_hbca_celltype)
ggsave(file.path(dir_res, "barplot.celltypes_TNBC_vs_HBCA.pdf"), p2, width = 6, height = 5, useDingbats = FALSE)

#------ Copykat vs HBCA ------
p3 <- ruok::qtable_scatter(copykat_res_df$copykat_pred_tirosh, rctd_res_df_hbca$first_type) +
    labs(x = "CopyKat (tirosh)", y = "HBCA") +
    scale_color_manual(values = pal_copykat)
ggsave(file.path(dir_res, "scatter.celltypes_Copykat_vs_HBCA.pdf"), p3, width = 3.5, height = 5, useDingbats = FALSE)
p4 <- ruok::qbarplot_table_catx(copykat_res_df$copykat_pred_tirosh, rctd_res_df_hbca$first_type, name_h = "CopyKat (tirosh)", name_v = "HBCA") +
    scale_fill_manual(values = pal_hbca_celltype)
ggsave(file.path(dir_res, "barplot.celltypes_Copykat_vs_HBCA.pdf"), p4, width = 3.5, height = 5, useDingbats = FALSE)

#------ Copykat vs HBCA ------
p3 <- ruok::qtable_scatter(copykat_res_df$copykat_pred_tirosh, rctd_res_df_hbca$first_type) +
    labs(x = "CopyKat (tirosh)", y = "HBCA") +
    scale_color_manual(values = pal_copykat)
ggsave(file.path(dir_res, "scatter.celltypes_Copykat_vs_HBCA.pdf"), p3, width = 3.5, height = 5, useDingbats = FALSE)
p4 <- ruok::qbarplot_table_catx(copykat_res_df$copykat_pred_tirosh, rctd_res_df_hbca$first_type, name_h = "CopyKat (tirosh)", name_v = "HBCA") +
    scale_fill_manual(values = pal_hbca_celltype)
ggsave(file.path(dir_res, "barplot.celltypes_Copykat_vs_HBCA.pdf"), p4, width = 3.5, height = 5, useDingbats = FALSE)

#------ Copykat vs TNBC ------
p5 <- ruok::qtable_scatter(copykat_res_df$copykat_pred_tirosh, rctd_res_df_tnbc$first_type) +
    labs(x = "CopyKat (tirosh)", y = "TNBC") +
    scale_color_manual(values = pal_copykat)
ggsave(file.path(dir_res, "scatter.celltypes_Copykat_vs_TNBC.pdf"), p5, width = 3.5, height = 5, useDingbats = FALSE)
p6 <- ruok::qbarplot_table_catx(copykat_res_df$copykat_pred_tirosh, rctd_res_df_tnbc$first_type, name_h = "CopyKat (tirosh)", name_v = "TNBC") +
    scale_fill_manual(values = pal_celltypes)
ggsave(file.path(dir_res, "barplot.celltypes_Copykat_vs_TNBC.pdf"), p6, width = 3.5, height = 5, useDingbats = FALSE)

#------ TNBC vs Copykat ------
p7 <- ruok::qtable_scatter(rctd_res_df_tnbc$first_type, copykat_res_df$copykat_pred_tirosh) +
    labs(x = "TNBC", y = "CopyKat (tirosh)") +
    scale_color_manual(values = pal_celltypes)
ggsave(file.path(dir_res, "scatter.celltypes_TNBC_vs_Copykat.pdf"), p7, width = 6, height = 5, useDingbats = FALSE)
p8 <- ruok::qbarplot_table_catx(rctd_res_df_tnbc$first_type, copykat_res_df$copykat_pred_tirosh, name_h = "TNBC", name_v = "CopyKat (tirosh)") +
    scale_fill_manual(values = pal_copykat)
ggsave(file.path(dir_res, "barplot.celltypes_TNBC_vs_Copykat.pdf"), p8, width = 6, height = 5, useDingbats = FALSE)

#------------------ ~~~ Finalize cell types ~~~ --------------------
cli_h1("Finalize cell types")

sr3 <- AddMetaData(sr3, metadata = rctd_res_df_hbca$first_type, col.name = "celltype_HBCA")
sr3 <- AddMetaData(sr3, metadata = rctd_res_df_tnbc$first_type, col.name = "celltype_TNBC")
sr3 <- AddMetaData(sr3, metadata = copykat_res_df$copykat_pred_tirosh, col.name = "copykat")

css_celltype <- sr3$celltype_TNBC %>% as.character()
# Refine tumor cells
if (!sum(sr3$copykat == "aneuploid") == 0) {
    ## Ideal case:
    ## Tumor=Tumor & Epithelial cells (the HBA lumsec, basal, lumhr) & Aenuploid cells (Copykat's aneuploid)
    idx <- sr3$celltype_TNBC == "Tumor" & c(c(!sr3$celltype_HBCA %in% c("basal", "lumhr", "lumsec")) | sr3$copykat != "aneuploid")
} else {
    ## Non-ideal case:
    ## Tumor=Tumor & Epithelial cells (the HBA lumsec, basal, lumhr)
    cli_alert_danger("Copykat aneuploid cells not found, using HBCA to refine tumor cells")
    idx <- sr3$celltype_TNBC == "Tumor" & c(c(!sr3$celltype_HBCA %in% c("basal", "lumhr", "lumsec")))
}
css_celltype[idx] <- "LOWCONF"

css_celltype <- factor(css_celltype, levels = names(pal_celltypes))
css_celltype <- forcats::fct_drop(css_celltype)
table(css_celltype, sr3$celltype_TNBC)
names(css_celltype) <- Cells(sr3)
css_celltype_df <- enframe(css_celltype, "Barcode", "celltypes")
head(css_celltype_df)
write_rds(css_celltype_df, file.path(dir_res, "celltype_final.rds"))
write.csv(css_celltype_df, file.path(dir_res, "celltype_final.csv"))

sr3 <- AddMetaData(sr3, metadata = css_celltype, col.name = "celltypes_css")

#------------------ ~~~ Spatial plot ~~~ --------------------
print(ncol(sr3))
sp_xy_ratio <- get_visium_xy_ratio(sr3)
sp_x_len <- GetTissueCoordinates(sr3)$x %>%
    range() %>%
    diff()
sp_y_len <- GetTissueCoordinates(sr3)$y %>%
    range() %>%
    diff()
print(c(sp_x_len, sp_y_len, sp_xy_ratio))
pt.size.factor <- suggest_seurat_visium_spatial_pt_size(sr3)
pdf_basic_inch <- 5
pdf_width <- pdf_basic_inch
pdf_height <- pdf_basic_inch
# pt.size.factor=5 will good at evenly filling the image with points
if (sp_xy_ratio > 1) {
    pdf_height <- pdf_height * sp_xy_ratio
} else {
    pdf_width <- pdf_width * sp_xy_ratio
}
print(c(pdf_width, pdf_height))

cli_h1("Spatial plot")

for (z in c("celltype_HBCA", "celltype_TNBC", "copykat", "celltypes_css")) {
    cli_alert_info(glue("Plot spatial dimplot for {z}"))
    pal_use <- switch(z,
        celltype_HBCA = pal_hbca_celltype,
        celltype_TNBC = pal_celltypes,
        copykat = pal_copykat,
        celltypes_css = pal_celltypes
    )
    p <- SpatialDimPlot(
        sr3,
        group.by = z, pt.size.factor = pt.size.factor,
        label = F, alpha = 1,
    ) + theme(aspect.ratio = sp_xy_ratio)
    p <- p + scale_fill_manual(values = pal_use, name = z)
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

    p <- p + scale_color_manual(values = pal_use, name = z)
    ggsave(
        filename = file.path(dir_res, sprintf("dimplot.%s.pdf", z)),
        plot = p,
        width = 4,
        height = 4, useDingbats = F
    )
}

if (T) {
    sr3$hi <- sr3$celltypes_css == "Tumor"
    p <- SpatialDimPlot(
        sr3,
        # cells.highlight = Cells(sr3)[which(sr3$hi)],
        # cols.highlight = c("red", '#bebebe00'),
        group.by = "hi",
        pt.size.factor = pt.size.factor,
        label = F,
    ) + theme(aspect.ratio = sp_xy_ratio)
    p <- p + scale_fill_manual(values = c("#bebebe0e", "red"), name = "Tumor cells")
    ggsave(
        filename = file.path(dir_res, "spatial.Tumor_cells.pdf"),
        plot = p + rremove("legend") + ggtitle(sprintf("%s Tumor cells (n=%d)", sample_name, sum(sr3$hi, na.rm = TRUE))),
        width = pdf_width + 3,
        height = pdf_height, useDingbats = F
    )
    p <- DimPlot(
        sr3,
        group.by = "hi",
        label = T, raster = T, pt.size = 3
    ) + theme(aspect.ratio = 1)
    p <- p + scale_color_manual(values = c("grey", "red"), name = "Tumor cells")
    ggsave(
        filename = file.path(dir_res, "dimplot.Tumor_cells.pdf"),
        plot = p + rremove("legend") + ggtitle(sprintf("%s Tumor cells (n=%d)", sample_name, sum(sr3$hi, na.rm = TRUE))),
        width = 4,
        height = 4, useDingbats = F
    )
    sr3$hi <- NULL
}

#------------------ ~~~ Copykat CNA heatmap ~~~ --------------------
# stop('As expected')
cli_h1("Copykat CNA heatmap")
source("/volumes//USR1/yyan/project/tumor_plasticity/core/GTA/package_general_cna_analysis.R")
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
# if (!file.exists(file.path(dir_res, "copykat_mix.sce.rds"))) {
if (TRUE) {
    f_copykat_gta <- file.path(
        dirname(f_in),
        "copykat_mix/copykat_mix.sce.rds"
    )
    gta_mix <- read_rds(f_copykat_gta)
    print(colnames(colData(gta_mix)))

    str(intersect(make.names(rownames(copykat_res_df)), colnames(gta_mix)))

    dim(gta_mix)
    sr3$in_gta <- make.names(Cells(sr3)) %in% colnames(gta_mix)
    table(sr3$in_gta, sr3$celltypes_css)


    if (!identical(colnames(gta_mix), Cells(sr3))) {
        cli::cli_alert_warning("[expected] Mismatch in cell names between CopyKat and Seurat object.")
        tmp <- intersect(make.names(Cells(sr3)), colnames(gta_mix))
        str(tmp)
        print(c("shared_cells=", length(tmp), "copykat_cells" = length(colnames(gta_mix)), "object_cells" = length(Cells(sr3))))
        if (length(tmp) / ncol(gta_mix) < 0.1) {
            stop("Less than 10% of cells in Seurat object are found in CopyKat results. Something must be wrong [error]")
        }
        gta_mix <- gta_mix[, tmp]
        rm(tmp)
    }

    idx <- match(colnames(gta_mix), make.names(Cells(sr3)))
    for (z in c("celltype_TNBC", "celltype_HBCA", "celltypes_css")) {
        tmp <- as.character(sr3@meta.data[[z]])[idx]
        if (z == "celltype_css") {
            tmp <- factor(tmp, levels = names(pal_celltypes))
        } else {
            tmp <- as.factor(tmp)
        }
        tmp <- forcats::fct_drop(tmp)
        colData(gta_mix)[, z] <- tmp
    }
    print(table(colData(gta_mix)$celltypes_css))

    write_rds(gta_mix, file.path(dir_res, "copykat_mix.sce.rds"))
} else {
    gta_mix <- read_rds(file.path(dir_res, "copykat_mix.sce.rds"))
}
#------ Viz ------
hm_obj_sc <- plot_heatmap_sc_manual(
    gta_mix,
    clip = c(-1, 1),
    cell_group_by = "celltypes_css",
    anno_rows_category = c(
        "celltype_TNBC",
        "celltype_HBCA",
        "copykat_pred_tirosh"
    ),
    row_anno_color_list = list(
        celltype_TNBC = pal_celltypes,
        celltype_HBCA = pal_hbca_celltype,
        copykat_pred_tirosh = pal_copykat
    )
)
pdf(file.path(dir_res, "copykat_heatmap3.celltypes_css.pdf"), width = 13, height = 10, useDingbats = FALSE)
draw(hm_obj_sc,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

cat("[DONE]")
