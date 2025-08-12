suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: clean and remove cells if they:
# - <= 10 nUMI
# - <= 10 nFeature
# - < 10 nportion.css_id
# - misfolded layered areas based on H&E imaging (this is probably because of cutting)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 20241014
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
    library(fs)
    library(stringr)
    library(tibble)
    library(Seurat)
    library(pbapply)
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")

    adhoc_load_HE_whitelist <- function(dir_in) {
        f_list <- Sys.glob(file.path(dir_in, "*stats.csv"))
        if (length(f_list) == 0) {
            warning("no white list file is detected.")
            return(NULL)
        }

        o <- lapply(f_list, function(f) {
            message(f, "...")
            df <- read_csv(f, comment = "#", show_col_types = FALSE)
            cname <- df[["Cell ID"]]
            message(c(length(cname), " cells ...\n"))
            return(cname)
        })

        adhoc_load_HE_blacklist <- adhoc_load_HE_whitelist
        o <- unique(as.character(unlist(o)))
        return(o)
    }
})
cmdargs <- commandArgs(trailingOnly = TRUE)
print(cmdargs)


if (length(cmdargs) > 0) {
    sample_name <- cmdargs[1]
} else {
    sample_name <- "ART266"
    sample_name <- "ART247"
}

#------ parameters ------

if (sample_name %in% c(
    "ART23", "ART10", "ART312", "ART3122", "ART311", "ART305", "ART304",
    "ART18", "ART31", "ART40", "ART43", "ART65"
)) {
    message("Xenium with OCT sapmles")
    param_nCount_Xenium_cutoff <- 10
    param_nFeature_Xenium_cutoff <- 10
} else {
    message("Xenium with FFPE samples")
    param_nCount_Xenium_cutoff <- 3
    param_nFeature_Xenium_cutoff <- 3
}

cli_rule(sample_name)

#------------------ ~~~ START ~~~ --------------------
cli_h1('START')
f_sp0 <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium/data/",
    sample_name,
    "nonbinarized_pca/xenium_nonbinarized_pca.seurat.rds"
)

f_css_id <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium/data",
    sample_name,
    "nonbinarized_pca",
    "TransferLabelEval.ref_RNA.query_Xenium.Consensus",
    "celltype", "Consensus_TransferData_css_id_full.dataframe.rds"
)

dir_he_whitelist <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium/HE_image",
    "xenium_HE_cell_whitelist", sample_name
)
dir_he_blacklist <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium/HE_image",
    "xenium_HE_cell_blacklist", sample_name
)

dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium/data/",
    sample_name, "cleaned1"
)
dir_create(dir_res)

dir_diagnose <- file.path(dir_res, "diagnose_cleaning")
fs::dir_create(dir_diagnose)
#------------------ ~~~ Read ~~~ --------------------
cli_h1('Read')
sp <- read_rds(f_sp0)
str(Cells(sp))

#------ basic metrics of nUMI and nFeatures ------
cnames_pass_basic_qc <- Cells(sp)[sp$nCount_Xenium > param_nCount_Xenium_cutoff & sp$nFeature_Xenium > param_nFeature_Xenium_cutoff]
str(cnames_pass_basic_qc)

#------------------ ~~~ Density removing outlier tissues ~~~ --------------------
## [optional]
## helps to determine the H&E whitelist
if (FALSE) {
    ## too slow
    df_coord <- GetTissueCoordinates(sp)
    head(df_coord)
    library(dbscan)
    library(parallelDist)
    df_coord <- tibble::column_to_rownames(df_coord, "cell")
    spatial_dist <- parallelDist(as.matrix(df_coord[, c("x", "y")]),
        method = "euclidean"
    )
    class(spatial_dist)
    set.seed(1026)
    dbscan_obj <- hdbscan(spatial_dist, minPts = 5)
}


#------ Consensus id (css_id_full) result ------
cli_h2('Consensus id (css_id_full) result')
df_css_id <- read_rds(f_css_id)
tmp <- match(Cells(sp), df_css_id$cellname)
df_css_id <- df_css_id[tmp, ]
df_css_id$cellname <- Cells(sp)
cells_to_css <- df_css_id$css_id_full[tmp]
is_bad_cell <- cells_to_css == "LOWCONF"
is_bad_cell[is.na(cells_to_css)] <- FALSE
# table(is_bad_cell)
cnames_pass_cssid <- df_css_id$cellname[!is_bad_cell]
str(cnames_pass_cssid)

print(table(rownames(df_css_id) %in% cnames_pass_basic_qc, df_css_id$css_id_full))

#------ HE blacklist ------
cli_h2('HE blacklist')
cnames_HE_blacklist <- adhoc_load_HE_whitelist(dir_he_blacklist)
if (is.null(cnames_HE_blacklist)) {
    message("no any black list so use all cells")
    cnames_HE_blacklist <- NULL
}
str(cnames_HE_blacklist)

#------ HE whitelist ------
cli_h2('HE whitelist')
cnames_HE_whitelist <- adhoc_load_HE_whitelist(dir_he_whitelist)
str(cnames_HE_whitelist)
if (is.null(cnames_HE_whitelist)) {
    message("no any white list so use all cells")
    cnames_HE_whitelist <- Cells(sp)
}

print(table(HE_white_list = Cells(sp) %in% cnames_HE_whitelist))
cnames_in_whitelist <- intersect(Cells(sp), cnames_HE_whitelist)
str(cnames_in_whitelist)

cnames_pass_HE <- intersect(setdiff(Cells(sp), cnames_HE_blacklist), cnames_in_whitelist)
str(cnames_pass_HE)
#------------------ ~~~ Report ~~~ --------------------
cli_h1('Report')
cnames_clean <- intersect(Cells(sp), cnames_pass_cssid) %>%
    intersect(., cnames_HE_whitelist) %>%
    intersect(., cnames_pass_basic_qc)
str(cnames_clean)
library(UpSetR)
p <- UpSetR::upset(UpSetR::fromList(
    list(
        raw = Cells(sp),
        pass_css = cnames_pass_cssid,
        pass_HE = cnames_pass_HE,
        pass_basic = cnames_pass_basic_qc,
        clean = cnames_clean
    )
))
pdf(file.path(dir_diagnose, "plot.upsetr.pdf"), width = 5, height = 5, useDingbats = F)
print(p)
dev.off()

report_df <- data.frame(
    raw        = length(Cells(sp)),
    pass_HE    = length(cnames_pass_HE),
    pass_basic = length(cnames_pass_basic_qc),
    pass_css   = length(cnames_pass_cssid),
    clean      = length(cnames_clean)
)
write_rds(report_df, file.path(dir_diagnose, "report.ncells_toreach_clean.rds"))
write_csv(report_df, file.path(dir_diagnose, "report.ncells_toreach_clean.csv"))


#------------------ ~~~ Attach the css id as celltype ~~~ --------------------
cli_h1('Attach the css id as celltype')

#------ gene dot plot before cleaning ------
cli_h2('gene dot plot before cleaning')
stopifnot(identical(df_css_id$cellname, Cells(sp)))
table(df_css_id$css_id)
sp <- AddMetaData(sp, metadata = df_css_id$css_id, col.name = "celltype_pre_qc")
if (T) {
    xmo <- sp
    ident_str <- "celltype_pre_qc"
    marker_suit <- "cellident"
    dir_snippet_viz <- dir_res
    f_df_marker <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/xenium/gene_probe.TME_celltypes/export.xenium.celltypes.csv"
    df_marker <- read.csv(f_df_marker)
    list_marker <- ruok::deframe_to_list(df_marker[, c("ident", "gene")])
    list_marker <- lapply(list_marker, function(xx) intersect(xx, rownames(xmo)))
    list_marker <- list_marker[intersect(levels(xmo@meta.data[, ident_str]), names(list_marker))]
    list_marker <- lapply(list_marker, sort)
    str(list_marker)
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet.dot_plot_genes.R")
}
sp$celltype_pre_qc <- NULL

#------ gene dot plot after cleaning ------
cli_h2('gene dot plot after cleaning')
stopifnot(identical(df_css_id$cellname, Cells(sp)))
sp <- AddMetaData(sp, metadata = df_css_id$css_id_full, col.name = "celltype")
sp <- subset(sp, cells = intersect(Cells(sp), cnames_clean))
str(Cells(sp))
library(forcats)
sp$celltype <- forcats::fct_drop(sp$celltype)

if (T) {
    xmo <- sp
    ident_str <- "celltype"
    marker_suit <- "cellident"
    dir_snippet_viz <- dir_res
    f_df_marker <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/xenium/gene_probe.TME_celltypes/export.xenium.celltypes.csv"
    df_marker <- read.csv(f_df_marker)
    list_marker <- ruok::deframe_to_list(df_marker[, c("ident", "gene")])
    list_marker <- lapply(list_marker, function(xx) intersect(xx, rownames(xmo)))
    list_marker <- list_marker[intersect(levels(xmo@meta.data[, ident_str]), names(list_marker))]
    list_marker <- lapply(list_marker, sort)
    str(list_marker)
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet.dot_plot_genes.R")
}

#------------------ ~~~ Export ~~~ --------------------
cli_h1('Export')
print(table(sp$celltype, useNA = "ifany"))
Idents(sp) <- "celltype"
message("export data")
write_xenium_rds(sp, dir_res, obj_type = "pass")

#------ viz idents ------
cli_h2('viz idents')
message("visulization")
dir_snippet_viz <- dir_res
xmo <- sp
pal_z <- pal_celltypes
viz_what <- "celltype"

ggsave(file.path(dir_res, sprintf("legend.%s.pdf", "celltype")),
    pal_to_ggplot(pal_z[levels(sp@meta.data[, "celltype"])], "celltype"),
    width = 3, height = 3, useDingbats = F
)

source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_categorical.R")


cli_alert_success('[DONE]')
#------ viz dotplot of cell type marker genes ------


if (F) {
    UMAPPlot(sp, group.by = "celltype")
    xxx <- Cells(sp)[Embeddings(sp, "umap")[, "UMAP_1"] < -8]
    xxx
    UMAPPlot(sp, cells.highlight = xxx)
    yyy <- GetTissueCoordinates(sp)
    head(yyy)
    table(yyy$hi)
    yyy$hi <- yyy$cell %in% xxx
    ggplot(yyy, aes(x = x, y = y)) +
        geom_point_rast(aes(color = hi), size = .1) +
        coord_equal() +
        scale_y_reverse() +
        scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "lightgrey"))
}
