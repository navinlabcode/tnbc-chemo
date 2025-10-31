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
lib_type <- NULL
if (length(cmdargs) > 0) {
    f_in <- cmdargs[1]
    sample_name <- cmdargs[2]
    lib_type <- cmdargs[3]
} else {
    ## Example Input
    sample_name <- "ART266"
    f_in <- glue("/volumes/USR1/yyan/project/tnbc_visium_hd/data/{sample_name}/ready.seurat.rds")
    lib_type <- "TNBC"
}

dir_res <- file.path(dirname(f_in), sprintf("celltype_infer_from_%s", lib_type))
fs::dir_create(dir_res)

#------ lib ------
f_lib_scRNA_HBCA <- "/core_users/users/ctang4/VisiumHD_Report/Annotation_Reference/Navin_HBCA_scRNA_10k_Ref_Ready.rds" ## Kindly provided by Chenling Tang
f_lib_scRNA_TNBC <- "/volumes/USR1/yyan/project/tnbc_visium_hd/lib/RCTD_Lib/RCTD_object.TNBC_celltypes.rds"

f_lib_scRNA <- switch(lib_type,
    "HBCA" = f_lib_scRNA_HBCA,
    "TNBC" = f_lib_scRNA_TNBC,
    stop("Unknown library type")
)
lib_scRNA <- read_rds(f_lib_scRNA)
print(class(lib_scRNA))

source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
pal_use <- switch(lib_type,
    "HBCA" = pal_hbca_celltype,
    "TNBC" = pal_celltypes
)
label_lvs <- names(pal_use)
print(label_lvs)
#------------------ ~~~ Readin ~~~ --------------------
cli_h1("Readin")
sr3 <- read_rds(f_in)
print(sr3)

f_rctd <- file.path(dir_res, "RCTD.rds")
if (file.exists(f_rctd)) {
    cli_h1("Readin RCTD")
    rctd_res <- read_rds(f_rctd)
} else {
    cli_h1("Run RCTD")
    # Run RCTD
    rctd_res <- Run_RCTD(sr3, lib_scRNA)
    saveRDS(rctd_res, file = f_rctd)
}

#------------------ ~~~ Viz ~~~ --------------------
cli_h1("Viz")
rctd_res_df <- rctd_res@results$results_df
print(colnames(rctd_res_df))
str(rctd_res_df)
f_rctd_sc <- file.path(dir_res, "RCTD.seurat_dataframe.rds")

try(all(rownames(rctd_res_df) == Cells(sr3)))

if (!identical(rownames(rctd_res_df), Cells(sr3))) {
    cli_alert_danger("The rownames of the RCTD results do not match the cell names in the Seurat object.")
    shared_ncells <- intersect(rownames(rctd_res_df), Cells(sr3))
    cli_alert_info(sprintf("%s / %s cells are shared between RCTD results and Seurat object.", 
        comma(length(shared_ncells)), comma(length(Cells(sr3)))))
    if (length(shared_ncells)/length(Cells(sr3)) < 0.1) {
        stop("Less than 10% of cells are shared between RCTD results and Seurat object; Something must be wrong.[error]")
    }
    idx <- match(Cells(sr3), rownames(rctd_res_df))
    rctd_res_df <- rctd_res_df[idx, , drop = FALSE]
    rownames(rctd_res_df) <- Cells(sr3)
    
    rctd_res_df$first_type[is.na(rctd_res_df$first_type)] <- "Unknown"
    rctd_res_df$scond_type[is.na(rctd_res_df$scond_type)] <- "Unknown"
    rctd_res_df$spot_class <- as.character(rctd_res_df$spot_class)
    rctd_res_df$spot_class[is.na(rctd_res_df$spot_class)] <- "Unknown"
}

rctd_res_df$first_type <- factor(rctd_res_df$first_type, levels = label_lvs)
rctd_res_df$scond_type <- factor(rctd_res_df$scond_type, levels = label_lvs)
rctd_res_df$first_type <- forcats::fct_drop(rctd_res_df$first_type)
rctd_res_df$scond_type <- forcats::fct_drop(rctd_res_df$scond_type)

write_rds(rctd_res_df, f_rctd_sc)
write.csv(rctd_res_df, str_replace_all(f_rctd_sc, "\\.rds$", ".csv"))

sr3 <- AddMetaData(sr3, metadata = rctd_res_df)


#------ Basic barplot of "spot_class"    "first_type"    "scond_type" ------
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
library(ruok)
library(patchwork)
for (z in c("spot_class", "first_type", "scond_type")) {
    if (!z %in% colnames(rctd_res_df)) {
        next()
    }
    cli_alert_info(glue("Plot {z}"))
    p <- ruok::qbarplot_table_cat(rctd_res_df[[z]], name_x = z)
    pp <- ruok::qbarplot_table_cat(rctd_res_df[[z]], name_x = z, do.prop.table = F)
    ## barplot with fill color
    pp <- ggplot(pp$data, aes(x = V, y = X, fill = X)) +
        geom_col() +
        geom_text(aes(label = V), position = position_stack(vjust = 0.9), size = 6 / .pt) +
        scale_y_discrete(limits = rev) +
        rremove("y.title") +
        labs(x = "num of cells")
    if (z %in% c("first_type", "scond_type")) {
        p <- p + scale_fill_manual(values = pal_use, name = z)
        pp <- pp + scale_fill_manual(values = pal_use, name = z)
    }

    ggsave(file.path(dir_res, glue("barplot.{z}.pdf")),
        patchwork::wrap_plots(p, pp, nrow = 1, guides = "collect", widths = c(1, 2)),
        width = 4, height = 2, useDingbats = F
    )
}

#------ first_type vs scond_type ------
p <- ruok::qtable_scatter(rctd_res_df$first_type, rctd_res_df$scond_type) +
    labs(x = "first_type", y = "scond_type") +
    scale_color_manual(values = pal_use)
ggsave(file.path(dir_res, "scatter.first_vs_second_type.pdf"), p, width = 6, height = 6, useDingbats = F)
#------ Dotplot of the known marker genes ------
if (lib_type == 'TNBC') {
    xmo <- sr3

    for (z in c("first_type", "scond_type")) {
        ident_str <- z

        marker_suit <- "celltypes"
        dir_snippet_viz <- dir_res
        f_df_marker <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/xenium/gene_probe.TME_celltypes/export.xenium.celltypes.csv"
        df_marker <- read.csv(f_df_marker)
        list_marker <- ruok::deframe_to_list(df_marker[, c("ident", "gene")])
        list_marker <- lapply(list_marker, function(xx) intersect(xx, rownames(xmo)))
        list_marker <- list_marker[intersect(levels(xmo[[]][, ident_str]), names(list_marker))]
        list_marker <- lapply(list_marker, sort)
        source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet.dot_plot_genes.R")
    }
    rm(xmo)
    rm(ident_str)
    rm(marker_suit)
    rm(df_marker)
    rm(list_marker)
}

#------ Spatial dimplot ------
sp_xy_ratio <- get_visium_xy_ratio(sr3)
pdf_basic_inch <- 5
pdf_width <- pdf_basic_inch
pdf_height <- pdf_basic_inch

if (sp_xy_ratio > 1) {
    pdf_height <- pdf_height * sp_xy_ratio
} else {
    pdf_width <- pdf_width * sp_xy_ratio
}
print(c(pdf_width, pdf_height))


for (z in c("spot_class", "first_type", "scond_type")) {
    if (!z %in% colnames(sr3[[]])) {
        next()
    }
    cli_alert_info(glue("Plot spatial dimplot for {z}"))
    p <- SpatialDimPlot(
        sr3,
        group.by = z, 
        # pt.size.factor = 4,
        label = F, alpha = 1,
    ) + theme(aspect.ratio = sp_xy_ratio)
    if (z %in% c("first_type", "scond_type")) {
        p <- p + scale_fill_manual(values = pal_use, name = z)
    }
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
    if (z %in% c("first_type", "scond_type")) {
        p <- p + scale_color_manual(values = pal_use, name = z)
    }
    ggsave(
        filename = file.path(dir_res, sprintf("dimplot.%s.pdf", z)),
        plot = p,
        width = 4,
        height = 4, useDingbats = F
    )
}


cat("[done]")
timestamp()
