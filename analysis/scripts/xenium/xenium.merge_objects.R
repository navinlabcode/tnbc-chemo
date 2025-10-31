# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge Xenium data objects
#
# Only use the Xenium assay
#
# Remove the spatial related assays (FOV) as they not needed
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(readr)
library(cli)
library(fs)
options <- commandArgs(trailingOnly = TRUE)

if (length(options) > 0) {
    dir_res <- options[1]
    f_in_opts <- options[2:length(options)]
} else {
    dir_res <- "/volumes/USR1/yyan/project/tnbc_xenium/_test"
    f_in_opts <- c(
        "/volumes/USR1/yyan/project/tnbc_xenium/data/ART40/nonbinarized_pca/celltyping.transferlabel.cca/objects_split_into_predicted.id/T/pass.sr3.rds",
        "/volumes/USR1/yyan/project/tnbc_xenium/data/ART43/nonbinarized_pca/celltyping.transferlabel.cca/objects_split_into_predicted.id/T/pass.sr3.rds",
        "/volumes/USR1/yyan/project/tnbc_xenium/data/ART65/nonbinarized_pca/celltyping.transferlabel.cca/objects_split_into_predicted.id/T/pass.sr3.rds"
    )
}

cli::cli_ol(f_in_opts)
dir_create(dir_res)

library(stringr)
sample_names_parsed <- stringr::str_extract(f_in_opts, "ART[0-9]+")
names(f_in_opts) <- sample_names_parsed

if (sum(!file_exists(f_in_opts)) > 0) {
    message("These input files did not exist:")
    cli::cli_ol(f_in_opts[which(!file_exists(f_in_opts))])
}

f_in_opts <- f_in_opts[file_exists(f_in_opts)]

if (length(f_in_opts) == 0) {
    stop("No solid inputs at all.")
}
if (length(f_in_opts) == 1) {
    stop("Only 1 input so no need to run merge.")
}

f_in_opts <- as.list(f_in_opts)

#------------------ ~~~ Rock ~~~ --------------------
rds_list <- lapply(names(f_in_opts), function(x) {
    o <- read_rds(f_in_opts[[x]])
    if ("FOV" %in% class(try(o[["fov"]]))) {
        o <- DietSeurat(o, assays = "Xenium")
        o[["fov"]] <- NULL
    }
    ## Add a `sample` column
    if (!"sample" %in% colnames(o@meta.data)) {
        o$sample <- x
    }
    ## Add the xenium panel information
    o$xenium_panel <- "x5k"
    if (x %in% paste0("ART", c(10, 23, 304, 305, 311, 312))) {
        o$xenium_panel <- "x5k"
    }
    if (x %in% paste0("ART", c(18, 31, 40, 43, 65))) {
        o$xenium_panel <- "x5kPlus"
    }

    return(o)
})
names(rds_list) <- names(f_in_opts)
print(rds_list)


str_col_use <- Reduce(intersect, lapply(rds_list, function(x) colnames(x@meta.data)))
# str_col_use <- c('major_groups', 'aneuploidy_tri_type',
#                  'sample_id', 'samplename',
#                  'patient_id', 'chemistry',
#                  'subject_id', 'patient',
#                  'RCB_status', 'PCR_status')
str_col_use <- str_col_use[!str_detect(str_col_use, "pred")]
str_col_use <- str_col_use[!str_detect(str_col_use, "snn_res")]
str_col_use <- setdiff(str_col_use, "seurat_clusters")
print(str_col_use)

rds_list <- lapply(rds_list, function(x) {
    df <- x@meta.data
    df <- df[, str_col_use, drop = F]
    x@meta.data <- df
    return(x)
})

rds_list <- lapply(rds_list, function(o) {
    o <- RenameCells(o, new.names = paste0(o$sample, "_", Cells(o)))
    o
})

cli_h1("Merge begins")
timestamp()
rds <- merge(rds_list[[1]],
    y = rds_list[2:length(rds_list)],
    project = "artemis_xenium"
)
timestamp()
print(rds)

cli_h1("Export data")
print(head(rds@meta.data))
write_rds(x = rds, path = file.path(dir_res, "pass.sr3.rds"))
write_rds(x = rds@meta.data, path = file.path(dir_res, "sr3_metadata.df.rds"))

timestamp()
cat("[done in R]")
