suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Create the raw Seurat object for Visium HD data
#       and save it as a RDS file.
# If the demultiplexing CSV file is provided, the code will
#       demultiplex the data and save the demultiplexed object.
# Seurat_5.1.0; R version 4.4.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025-05-05
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
    library(hdf5r)
    library(arrow)
})
cmdargs <- commandArgs(trailingOnly = TRUE)

vs_binsize <- 0
f_demultiplex_csv <- NULL
FIXED_VS_BINSIZE_INPUT <- 32 # By default, this folder always has the 8um and 16um data
if (length(cmdargs) > 0) {
    sample_name <- cmdargs[1]
    vs_binsize <- as.numeric(cmdargs[2])
    print(cmdargs)
} else {
    sample_name <- "ART223"
    vs_binsize <- 32

}

vs_dir_in <- file.path(
    "/volumes/USR1/yyan/project/tnbc_visium_hd/data0",
    sample_name,
    sprintf("binsize_%s", vs_binsize),
    "loupe_outs"
)


cli_alert_info("Sample name: {sample_name}")
cli_alert_info("Bin size: {vs_binsize}")
cli_alert_info("Input directory: {vs_dir_in}")

dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_visium_hd/data0",
    sample_name,
    sprintf("binsize_%s", vs_binsize)
)
fs::dir_create(dir_res)


#------------------ ~~~ Read-in ~~~ --------------------
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
cli_h1("Read-in")

f_obj <- file.path(dir_res, "pass.seurat.rds")
if (!file.exists(f_obj)) {
    cli_alert_info("Read-in the Visium HD data")
    object <- Load10X_Spatial(data.dir = vs_dir_in, bin.size = c(vs_binsize))
    print(object)
    orig_assay_name <- DefaultAssay(object)
    ## rename the assay to "Spatial"
    object <- RenameAssays(object, orig_assay_name, "Spatial")

    write_seurat(object = object, dir_out = dir_res, obj_type = "pass")
} else {
    cli_alert_info("The Seurat object is already created; so load it")
    object <- read_rds(f_obj)
}

p <- SpatialPlot(object, features = "nCount_Spatial")
ggsave(
    filename = file.path(dir_res, "nCount_Spatial.pdf"),
    plot = p,
    width = 7,
    height = 7, useDingbats = F
)

f_umi <- file.path(dir_res, "umi_count.matrix.rds")
if (!file.exists(f_umi)) {
    cli_alert_info("Save the UMI count matrix")
    object %>%
        GetAssayData(slot = "counts") %>%
        readr::write_rds(f_umi)
} else {
    cli_alert_info("The UMI count matrix is already created; so load it")
    # object@assays$Spatial@counts <- read_rds(f_umi)
}

cat('DONE\n')
timestamp()
