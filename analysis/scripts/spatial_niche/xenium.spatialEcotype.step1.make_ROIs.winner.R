suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Making ROIs with cell states frequencies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-12-17
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
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.spatial_ecotypes.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
    library(InSituCor)
    library(MatrixGenerics)
    library(pbapply)
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    obj_suit <- cmdargs[1]
    radius_final_use <- as.numeric(cmdargs[2])

} else {
    obj_suit <- "ALL"
    radius_final_use <- 30 ## !!!parameter!!!
}


#------------------ ~~~ Pre-requisites ~~~ --------------------
cli_h1("Pre-requisites")


if (T) {
    samples_to_collect <- c(
        "ART10", "ART23", "ART312", "ART3122", "ART311", "ART305", "ART304",
        "ART18", "ART31", "ART40", "ART43", "ART65",
        "ART94", "ART133",
        "ART117", "ART217", "ART258", "ART272", "ART232", "ART282", "ART289", "ART122", "ART250", "ART92",
        "ART153", "ART194", "ART238", "ART257", "ART104", "ART279", "ART30", "ART71", "ART247", "ART271",
        "ART155", "ART170", "ART219", "ART236", "ART277", "ART202", "ART223", "ART235", "ART266", "ART276"
    )
    print(length(samples_to_collect))
    sample_cohort_name <- "data_merged_N44"
}
#------------------ ~~~ Inputs ~~~ --------------------
cli_h1("Inputs")

#------ Load data ------
dir_proj <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium", 
    sample_cohort_name, 
    "spatial_ecotype_winner",
    obj_suit
)
df_meta <- read_rds(file.path(dir_proj, "inputs", "dataframe.cellmeta.rds"))
df_value <- read_rds(file.path(dir_proj, "inputs", "dataframe.value.rds"))

#------------------ ~~~ Output ~~~ --------------------
cli_h1("Output")
# dir_res <- file.path(dir_proj, "spEcotypes")
dir_res <- file.path(dir_proj, "scNicheRadius")

fs::dir_create(dir_res)

#------ colors ------
pal_sample <- init_pal_d(samples_to_collect, "bear")
# show_col(pal_sample)
ggsave(file.path(dir_res, "pal_sample.pdf"), pal_to_ggplot(pal_sample, "sample"), width = 1.5, height = 3, useDingbats = F)


#------------------ ~~~ START ~~~ --------------------
cli_h1("START")
# !!! Inputs !!!
# - df_meta
# - df_value

stopifnot(identical(rownames(df_meta), rownames(df_value)))


if ("cellname" %in% colnames(df_value)) {
    df_value$cellname <- NULL
} # remove cellname column
mat_value <- as.matrix(df_value)
if(is.na(max(mat_value))) {
    warning("NA in the indicator matrix")
}
rm(mat_value)

feature_opts <- colnames(df_value)
cli_ol(feature_opts) # 60 features


#------------------ ~~~ Spatial neighborhood matrix ~~~ --------------------
cli_h1("Spatial neighborhood matrix")
# n_rand_cells_per_sample <- 1e5 ## !!!parameter!!! [DEPRECATED; downsampling is wrong]
print(table(df_meta$sample))

df_meta$sample_index <- as.integer(factor(df_meta$sample))

#------ Check: distance range per samples ------
if (T) {
    cli_h3("Check: distance range per samples")
    df_xy_range <- df_meta %>%
        group_by(sample) %>%
        summarize(
            min_x = min(coord_x_adj),
            max_x = max(coord_x_adj),
            min_y = min(coord_y_adj),
            max_y = max(coord_y_adj)
        ) %>%
        mutate(
            x_range = max_x - min_x,
            y_range = max_y - min_y
        )
    print(df_xy_range)
    write_csv(df_xy_range, file.path(dir_res, "check.sample_xy_range.csv"))
}

#------ Justify Radius ------
## Understand how choise of radius affects the number of cells per ROI
radius_to_search <- c(1, 5, 10, 25, 30, 50, 75, 100, 200)
if (T) {
    ## build a list of spatial NB matrix depending on the radius
    for (radius_use in radius_to_search) {
        cli_h3(c("radius_use = ", radius_use))
        f_sp_nb_mat <- file.path(dir_res, sprintf("spatial_nb_mat.R%s.rds", radius_use))
        if (file.exists(f_sp_nb_mat)) {
            cat("Loading spatial NB matrix...")
            sp_nb_mat <- read_rds(f_sp_nb_mat)
        } else {
            cat("Calculating spatial NB matrix...")
            sp_nb_mat <- radiusBasedGraph(
                x = df_meta$coord_x_adj,
                y = df_meta$coord_y_adj,
                R = radius_use,
                subset = df_meta$sample_index
            )
            write_rds(sp_nb_mat, f_sp_nb_mat)
        }
        cat("[done]\n")
    }
}
print(dim(sp_nb_mat))
write_rds(rownames(df_meta), file.path(dir_res, "sp_nb_mat.rowcolnames.rds"))

sp_nb_mat_rowcolnames <- rownames(df_meta)
str(sp_nb_mat_rowcolnames)
for (radius_use in radius_to_search) {
    cli_h3(c("radius_use = ", radius_use))

    ## calculate ncells of ROIs per sample depending on the radius
    f_sp_nb_mat <- file.path(dir_res, sprintf("spatial_nb_mat.R%s.rds", radius_use))
    sp_nb_mat <- read_rds(f_sp_nb_mat)

    n_cells_per_roi <- MatrixGenerics::rowSums(sp_nb_mat > 0)
    names(n_cells_per_roi) <- sp_nb_mat_rowcolnames
    print(quantile(n_cells_per_roi))
    df_ncells_per_roi_persample <- data.frame(
        sample = df_meta[names(n_cells_per_roi), "sample"],
        n_cells_per_roi = n_cells_per_roi
    )
    df_ncells_of_roi_persample_summary <- df_ncells_per_roi_persample %>%
        group_by(sample) %>%
        summarize(n_cells = round(mean(n_cells_per_roi)))
    write_rds(
        df_ncells_per_roi_persample,
        file.path(dir_res, sprintf("check.R%s.ncells_per_roi_persample.rds", radius_use))
    )
    write_rds(
        df_ncells_of_roi_persample_summary,
        file.path(dir_res, sprintf("check.R%s.ncells_of_roi_persample_summary.rds", radius_use))
    )
}

if (T) {
    ## plot average nCells of ROIs per sample depending on the radius
    df_ncells_of_roi_persample_summary_combo <- lapply(radius_to_search, function(radius_use) {
        read_rds(
            file.path(dir_res, sprintf("check.R%s.ncells_of_roi_persample_summary.rds", radius_use))
        ) %>% mutate(radius_use = radius_use)
    }) %>% bind_rows()
    p <- ggplot(data = df_ncells_of_roi_persample_summary_combo, aes(x = radius_use, y = n_cells)) +
        geom_point(size = 2, alpha = .5, aes(color = sample)) +
        geom_line(aes(color = sample)) +
        scale_color_manual(values = pal_sample)
    # p <- p + rremove("legend") + scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl")
    p <- p + rremove("legend")
    p <- p + scale_x_log10() + annotation_logticks(sides = "b")
    log10p <- function(x) {
        log10(x + 1)
    }
    log10p_inv <- function(x) {
        10^x - 1
    }
    scale_trans_log10p <- scales::trans_new("log10p", transform = log10p, inverse = log10p_inv)
    p <- p + scale_y_continuous(transform = scale_trans_log10p) + annotation_logticks(sides = "l")

    df_tmp_label <- df_ncells_of_roi_persample_summary_combo %>%
        group_by(radius_use) %>%
        summarize(n_cells = round(mean(n_cells))) %>%
        mutate(label = glue("{n_cells}"))
    p <- p + geom_text(data = df_tmp_label, aes(x = radius_use, y = n_cells, label = label))
    ggsave(
        file.path(dir_res, "check_plot.radius_vs_nCells_of_ROIs_per_sample.pdf"),
        p,
        width = 5, height = 3,
        useDingbats = F
    )
}


#------------------ ~~~ Make ROIs ~~~ --------------------
cli_h1("Make ROIs")
### pick the proper radius
radius_use <- radius_final_use
cli_alert_info(c("radius determined = ", radius_use))
### load the spatial NB matrix
sp_nb_mat <- read_rds(file.path(dir_res, sprintf("spatial_nb_mat.R%s.rds", radius_use)))
sp_nb_mat_rowcolnames <- read_rds(file.path(dir_res, "sp_nb_mat.rowcolnames.rds"))

dir_res <- file.path(dir_res, sprintf("R%s", radius_use))
fs::dir_create(dir_res)

dim(sp_nb_mat)

n_cells_per_roi <- MatrixGenerics::rowSums(sp_nb_mat > 0)
print(quantile(n_cells_per_roi))

print(table(df_meta$sample))
table(df_meta$sample) %>% enframe('sample', 'n_cells') %>% write.csv(file.path(dir_res, "check.sample_ncells.csv"), row.names = F)
table(df_meta$sample) %>% quick_summary_stats() %>% write.csv(file.path(dir_res, "check.sample_ncells.quick_stats_summary.csv"), row.names = F)
n_roi_per_sample <- round(median(table(df_meta$sample)))
print(n_roi_per_sample) #37,266

if (T) {
    ## full data
    df_ind <- as.matrix(df_value)
    stopifnot(identical(rownames(df_ind), sp_nb_mat_rowcolnames))
    sp_nb_mat_indicator <- (sp_nb_mat > 0) * 1
    full_composition_roi_ncell <- sp_nb_mat_indicator %*% df_ind
    n_cells_per_roi_full <- n_cells_per_roi
    full_composition_roi_freq <- full_composition_roi_ncell / n_cells_per_roi_full
    rownames(full_composition_roi_ncell) <- sp_nb_mat_rowcolnames
    rownames(full_composition_roi_freq) <- sp_nb_mat_rowcolnames
    write_rds(full_composition_roi_freq, file.path(dir_res, "result.full.composition_roi_freq.mat.rds"))
    write_rds(full_composition_roi_ncell, file.path(dir_res, "result.full.composition_roi_ncell.mat.rds"))
}

#------------------ ~~~ Make ROIs (downsampling) ~~~ --------------------
## Down-sampled data [do no use in the final analysis]
## Downsampling cells to make sure that
## 1) each sample is well represented,
## 2) both loose (low-ncell) and dense (high-ncells) ROIs are well represented.
n_cells_per_roi <- MatrixGenerics::rowSums(sp_nb_mat > 0)
names(n_cells_per_roi) <- sp_nb_mat_rowcolnames
# 1,785,903 cells in total
df_roi_info <- tibble(
    roi_name = sp_nb_mat_rowcolnames,
    n_cells = n_cells_per_roi,
    sample = df_meta[sp_nb_mat_rowcolnames, "sample"]
)
print(head(df_roi_info))
df_roi_info$n_cells %>% quick_summary_stats() %>% write.csv(file.path(dir_res, "check.full_roi_ncells.quick_stats_summary.csv"), row.names = F)

## randomly select n_roi_per_sample ROIs per sample
set.seed(1026)
df_roi_use <- df_roi_info %>%
    group_by(sample) %>%
    slice_sample(n = n_roi_per_sample) %>%
    ungroup()
table(df_roi_use$sample)

## todo: check the distribution of n_cells per ROI

roi_use <- df_roi_use$roi_name
str(roi_use)
roi_use <- intersect(sp_nb_mat_rowcolnames, roi_use)
write_rds(roi_use, file.path(dir_res, "roi_use.rds"))

## --- Count feature freq per ROI

calc_cell_identity_freq_per_ROI <- function(cellnames, df_ind, cell_identities_wanted) {
    if (!length(cellnames) > 1) {
        warning("Require > 1 cells")
        return(NULL)
    }
    if (!all(cellnames %in% rownames(df_ind))) {
        warning("Some cell names are not in the data frame")
        return(NULL)
    }
    if (missingArg(cell_identities_wanted)) {
        cell_identities_wanted <- colnames(df_ind)
    }
    cell_identities_wanted <- intersect(cell_identities_wanted, colnames(df_ind))

    df_ind <- df_ind[cellnames, cell_identities_wanted, drop = F]
    n_cells_this_roi <- length(cellnames)
    n_cells_each_feature <- colSums(df_ind, na.rm = T)
    freq_each_feature <- n_cells_each_feature / n_cells_this_roi
    return(freq_each_feature)
}

roi_use_idx <- match(roi_use, sp_nb_mat_rowcolnames)

if (T){
    ## much faster
    ## composition_roi_freq = NB_graph * df_indicator
    sp_nb_mat_use <- sp_nb_mat[roi_use_idx, roi_use_idx, drop = F]
    sp_nb_mat_use <- sp_nb_mat_use > 0
    sp_nb_mat_use <- sp_nb_mat_use * 1
    df_ind_use <- df_value[roi_use, ] %>% as.matrix()
    composition_roi_ncell <- sp_nb_mat_use %*% df_ind_use
    rownames(composition_roi_ncell) <- roi_use
    n_cells_per_roi_use <- n_cells_per_roi[roi_use]
    composition_roi_freq <- composition_roi_ncell / n_cells_per_roi_use
}

## some ROIs have 1 cell only (i.e., zero NB cells)
roi_avail <- rownames(composition_roi_freq)
stopifnot(all(rownames(composition_roi_freq) == names(composition_roi_ncell)))

if (T) {
    ## export results
    write_rds(composition_roi_freq, file.path(dir_res, "result.composition_roi_freq.mat.rds"))
    write_rds(composition_roi_ncell, file.path(dir_res, "result.composition_roi_ncell.mat.rds"))

    write_rds(roi_avail, file.path(dir_res, "result.roi_avail.rds"))
}

## plot ncells distribution of rois
if (T) {
    print(quantile(n_cells_per_roi_use))
    p <- ggplot(data = enframe(n_cells_per_roi_use), aes(x = value)) +
        geom_histogram(binwidth = 5, fill = "skyblue") +
        labs(x = "Number of cells per ROI", y = "Number of ROIs")
    ggsave(
        file.path(dir_res, "result_plot.ncells_distribution_of_rois.pdf"),
        p,
        width = 3, height = 2,
        useDingbats = F
    )
}
#------------------ ~~~ END ~~~ --------------------


cat("[done]")
timestamp()
