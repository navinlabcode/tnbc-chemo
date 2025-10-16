suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Creating inputs to run spatial ecotype analysis.
#
# Targeted Results:
# - a cells x features data frame
#   - Cancer MPs indicator: 1=yes, 0=no
#   - TME states indicator: 1=yes, 0=no
# - a cell meta.data data.frame: XY coordinates, sample, etc.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-12-17
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
suppressPackageStartupMessages({
    library(gtools)
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
    library(tidyr)
    library(ggrastr)
    library(scattermore)
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    library(Seurat)
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {} else {}

#------------------ ~~~ Inputs & Outputs ~~~ --------------------
cli_h1("Inputs & Outputs")
ident_is_what <- "cell_state_paper" ## FIXED

#------ Sample choices ------
if (TRUE) {
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

dir_proj <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium",
    sample_cohort_name,
    "spatial_ecotype_winner"
)
fs::dir_create(dir_proj)

#------------------ ~~~ START ~~~ --------------------
cli_h1("START")
col_wanted <- c("celltype", "sample", "cell_state_paper", "coord_x", "coord_y", "cellname")
col_wanted <- c(col_wanted, c("nCount_Xenium", "nFeature_Xenium")) ## technical factors
cancer_MP_method_wanted <- c("module") ## choices: module aucell ucell AMSREV
#------------------ ~~~ Load TME ~~~ --------------------
cli_h1("Load TME")
f_df_opts <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium",
    "data",
    samples_to_collect,
    "cleaned1",
    "cleaned2_TME",
    "pass.seurat.rds"
)
stopifnot(all(file.exists(f_df_opts)))
names(f_df_opts) <- samples_to_collect


f_res_TME <- file.path(dir_proj, "dataframe.TME.rds")
df_list <- lapply(samples_to_collect, function(s) {
    cat(which(samples_to_collect %in% s), "/", length(samples_to_collect), s, ": reading... ")
    obj_s <- read_rds(f_df_opts[s])

    cat("fetching meta.data ... ")
    df_s <- obj_s@meta.data
    if (!"cellname" %in% colnames(df_s)) {
        df_s$cellname <- Cells(obj_s)
    }
    if (!"sample" %in% colnames(df_s)) {
        df_s$sample <- s
    }

    cat("fetching spatial coordinates... ")
    df_xy_s <- GetTissueCoordinates(obj_s) # colnames are fixed: x y cell
    df_xy_s <- column_to_rownames(df_xy_s, "cell")
    if (!all(rownames(df_s) %in% rownames(df_xy_s))) {
        cat("WARNING: not all cells have spatial coordinates... ")
    } else {
        df_xy_s <- df_xy_s[rownames(df_s), ]
    }
    colnames(df_xy_s) <- paste0("coord_", colnames(df_xy_s))

    df_s <- cbind(df_s, df_xy_s)

    cat("adding sample prefix ... ")
    df_s$cellname <- paste0(as.character(df_s$sample), "_", df_s$cellname)
    rownames(df_s) <- df_s$cellname
    df_s <- df_s[, intersect(colnames(df_s), col_wanted)]

    cat("[done]\n")
    return(df_s)
})

df_TME <- do.call(rbind, df_list)

print(tail(df_TME, 2))
#                  nCount_Xenium nFeature_Xenium celltype         cellname sample
# ART65_ojdaboio-1           134             114        T ART65_ojdaboio-1  ART65
# ART65_ojkbhidb-1           236             190      Mye ART65_ojkbhidb-1  ART65
#                  cell_state_paper  coord_x  coord_y
# ART65_ojdaboio-1      NK-CD16high 3507.864 651.3602
# ART65_ojkbhidb-1             mast 3590.904 518.8488

write_rds(df_TME, f_res_TME)
write.csv(df_TME, paste0(f_res_TME, ".csv"))
nanoparquet::write_parquet(df_TME, paste0(f_res_TME, ".parquet"))

rm(df_list)
rm(f_df_opts)
#------------------ ~~~ Load Cancer Cells ~~~ --------------------
cli_h1("Load Cancer cells")
#------ Fetch meta.data and spatial coord ------
f_df_opts <- file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium",
    "data",
    samples_to_collect,
    "cleaned1", "objects_split_into_celltype", "Tumor", "pass.xenium.seurat.rds"
)
names(f_df_opts) <- samples_to_collect
stopifnot(all(file.exists(f_df_opts)))
# print(sum(file.exists(f_df_opts)))
# print(names(f_df_opts)[!file.exists(f_df_opts)])

f_res_cancer_meta <- file.path(dir_proj, "dataframe.Cancer_meta.rds")
if (!file_exists(f_res_cancer_meta)) {
    df_list <- lapply(samples_to_collect, function(s) {
        cat(which(samples_to_collect %in% s), "/", length(samples_to_collect), " ", s, ": reading... ")
        obj_s <- read_rds(f_df_opts[s])

        cat("fetching meta.data and spatial coordinates... ")
        df_s <- obj_s@meta.data
        if (!"cellname" %in% colnames(df_s)) {
            df_s$cellname <- Cells(obj_s)
        }
        if (!"sample" %in% colnames(df_s)) {
            df_s$sample <- s
        }

        df_xy_s <- GetTissueCoordinates(obj_s) # colnames are fixed: x y cell
        df_xy_s <- column_to_rownames(df_xy_s, "cell")
        if (!all(rownames(df_s) %in% rownames(df_xy_s))) {
            cat("WARNING: not all cells have spatial coordinates... ")
        } else {
            df_xy_s <- df_xy_s[rownames(df_s), ]
        }
        colnames(df_xy_s) <- paste0("coord_", colnames(df_xy_s))
        df_s <- cbind(df_s, df_xy_s)

        cat("adding sample prefix ... ")
        df_s$cellname <- paste0(as.character(df_s$sample), "_", df_s$cellname)
        rownames(df_s) <- df_s$cellname
        df_s <- df_s[, intersect(colnames(df_s), col_wanted)]

        cat("[done]\n")
        return(df_s)
    })

    df_cancer_meta <- do.call(rbind, df_list)

    head(df_cancer_meta, 2)
    tail(df_cancer_meta, 2)


    write_rds(df_cancer_meta, f_res_cancer_meta)
    write.csv(df_cancer_meta, paste0(f_res_cancer_meta, ".csv"))
    nanoparquet::write_parquet(df_cancer_meta, paste0(f_res_cancer_meta, ".parquet"))

    rm(df_list)
    rm(f_df_opts)
} else {
    df_cancer_meta <- read_rds(f_res_cancer_meta)
}

#------ Fetch MPs assignment ------
f_res_mp <- file.path(dir_proj, sprintf("dataframe.MP.%s.rds", cancer_MP_method_wanted))

df_mp <- read_rds(file.path(
    "/volumes/USR1/yyan/project/tnbc_xenium",
    sample_cohort_name,
    "celltype_Tumor",
    "AddModuleScore_MP",
    sprintf("MP_assignment_%sscore_cutoff_0.1_winner", cancer_MP_method_wanted),
    "data_logic.mat.rds"
))
head(df_mp)
stopifnot(all(rownames(df_mp) %in% rownames(df_cancer_meta)))
stopifnot(nrow(df_mp) == nrow(df_cancer_meta))

if (!identical(rownames(df_mp), rownames(df_cancer_meta))) {
    cli_alert_info("Polishing row order")
    idx <- match(rownames(df_cancer_meta), rownames(df_mp))
    df_mp <- df_mp[idx, ]
    rownames(df_mp) <- rownames(df_cancer_meta)
}

write_rds(df_mp, f_res_mp)
write.csv(df_mp, paste0(f_res_mp, ".csv"))
nanoparquet::write_parquet(df_mp, paste0(f_res_mp, ".parquet"))

# df_mp_score <- read_rds(
#     file.path(
#         "/volumes/USR1/yyan/project/tnbc_xenium",
#         sample_cohort_name,
#         "celltype_Tumor",
#         "AddModuleScore_MP",
#         "modulescore.df.rds"
#     )
# )
# print(colnames(df_mp_score))
# print(tail(df_mp_score, 2))


#------------------ ~~~ Combine TME and Cancer Cells ~~~ --------------------
cli_h1("Combine TME and Cancer MPs")

meta_colnames <- intersect(colnames(df_cancer_meta), colnames(df_TME))
print(meta_colnames)

#------ cancer meta and value ------
df_cancer_meta <- df_cancer_meta[, meta_colnames]

stopifnot(identical(rownames(df_mp), rownames(df_cancer_meta)))
head(df_mp, 2)
head(df_cancer_meta, 2)

#------ TME meta and value ------
df_TME_meta <- df_TME[, meta_colnames]

tme_ident_lvs <- levels(df_TME$cell_state_paper)
print(tme_ident_lvs)
print(table(df_TME$cell_state_paper))
df_TME_ind <- df_TME[, c("cellname", "cell_state_paper")]
df_TME_ind$ind <- 1
df_TME_ind <- pivot_wider(df_TME_ind, id_cols = "cellname", names_from = "cell_state_paper", values_from = "ind", values_fill = 0)
df_TME_ind <- df_TME_ind[, c("cellname", intersect(tme_ident_lvs, colnames(df_TME_ind)))]
# rownames(df_TME_ind) <- df_TME_ind$cellname
df_TME_ind <- df_TME_ind %>% column_to_rownames("cellname")

if (!identical(rownames(df_TME_meta), rownames(df_TME_ind))) {
    cli_alert_info("Poloshing row order")
    idx <- match(rownames(df_TME_meta), rownames(df_TME_ind))
    df_TME_ind <- df_TME_ind[idx, ]
    rownames(df_TME_ind) <- rownames(df_TME_meta)
}

#------ combine meta ------
df_meta <- rbind(df_cancer_meta, df_TME_meta)
#------ combine value ------

stopifnot(length(intersect(df_mp$cellname, rownames(df_TME_ind))) == 0)

# df_value <- full_join(df_mp, df_TME_ind, by='cellname') ## do not recommend

cancer_val_colnames <- setdiff(colnames(df_mp), "cellname")
tme_val_colnames <- setdiff(colnames(df_TME_ind), "cellname")

for (tmp in tme_val_colnames) {
    if (!tmp %in% colnames(df_mp)) {
        df_mp[[tmp]] <- 0
    }
}
for (tmp in cancer_val_colnames) {
    if (!tmp %in% colnames(df_TME_ind)) {
        df_TME_ind[[tmp]] <- 0
    }
}
rm(tmp)

stopifnot(identical(sort(colnames(df_mp)), sort(colnames(df_TME_ind))))
val_colnames <- colnames(df_mp)

df_value <- rbind(df_mp[, val_colnames], df_TME_ind[, val_colnames])
stopifnot(all(str_detect(rownames(df_value), "^ART")))
df_value$cellname <- rownames(df_value)

head(df_value, 2)
tail(df_value, 2)

#------------------ ~~~ Adjust spatial coordinates ~~~ --------------------
cli_h1("Adjust spatial coordinates")
n_samples <- length(unique(df_meta$sample))
pal_sample <- init_pal_d(samples_to_collect, "bear")
ggsave(file.path(dir_proj, "legend.sample.pdf"), pal_to_ggplot(pal_sample), width = 3, height = 3)

#------ original spatial coords ------

df_xy_range <- df_meta %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
        xmin = min(coord_x), xmax = max(coord_x),
        ymin = min(coord_y), ymax = max(coord_y)
    )
df_xy_range <- column_to_rownames(df_xy_range, "sample")
df_xy_range

#------ very big hjust ------
## avoid samples staying nearby. In reality, their distance should be infinite.
hjust_unit <- 3 * max(c(with(df_xy_range, xmax - xmin), with(df_xy_range, ymax - ymin)))
hjust_unit

#------ first, adjust each tissue to the origin ------
range(df_meta$coord_x)
range(df_meta$coord_y)
## x - xmin; y - ymin
df_meta$coord_x_adj_to_origin <- df_meta$coord_x - df_xy_range[df_meta$sample, "xmin"]
df_meta$coord_y_adj_to_origin <- df_meta$coord_y - df_xy_range[df_meta$sample, "ymin"]


#------ then, horitonally spread the samples ------
## x + hjust_unit * (i-1), where i is the sample index
## y remains unchanged
samples_to_collect
sample_index <- structure(seq_along(samples_to_collect), names = samples_to_collect)
df_meta$coord_x_adj <- df_meta$coord_x_adj_to_origin + hjust_unit * (sample_index[df_meta$sample] - 1)
df_meta$coord_y_adj <- df_meta$coord_y_adj_to_origin


#------------------ ~~~ Spatial distance matrix ~~~ --------------------
library(Rfast)
# Computer Memory Issues
# spdist_mat <- df_meta %>%
#     dplyr::select(coord_x_adj, coord_y_adj) %>%
#     Rfast::Dist(method = "euclidean")
head(df_meta)
stopifnot("sample" %in% colnames(df_meta))

## Run Dist per sample and then combine
## to-do


#------------------ ~~~ Export ~~~ --------------------
cli_h1("Export")
head(rownames(df_meta))
head(rownames(df_value))
stopifnot(identical(rownames(df_meta), rownames(df_value)))

write_rds(df_meta, file.path(dir_proj, "dataframe.cellmeta.rds"))
write_rds(df_value, file.path(dir_proj, "dataframe.value.rds"))


#------------------ ~~~ Demultiplex into different object_suit ~~~ --------------------
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")

df_meta <- read_rds(file.path(dir_proj, "dataframe.cellmeta.rds"))
df_value <- read_rds(file.path(dir_proj, "dataframe.value.rds"))
head(df_meta)
colnames(df_value)
pal_sample <- init_pal_d(df_meta$sample, "parade")


cancer_feature_names <- grep("M\\d{2}$", colnames(df_value), value = TRUE)
tme_feature_names <- setdiff(intersect(colnames(df_value), std_idents_cellstates), cancer_feature_names)

stopifnot(setdiff(colnames(df_value), c(cancer_feature_names, tme_feature_names)) == "cellname")


colnames(df_meta)
print(table(df_meta$celltype, useNA = "ifany"))
## Features may have NA values in some cells
sapply(df_value, function(x) sum(is.na(x))) %>%
    enframe() %>%
    filter(value > 0)
# 1 module_score_M10 86395
# 2 module_score_M11 86395

object_suit_opts <- c(
    # "TME_only",
    # paste0("TME_and_", cancer_feature_names),
    # "Tumor_only",
    # "TME_and_tumor",
    "ALL"
)
# cli_ol(object_suit_opts)

stopifnot(identical(df_meta$cellname, rownames(df_value)))
stopifnot(identical(rownames(df_meta), rownames(df_value)))

#------ lib ------
if (T) {
    lib_patient_info <- read_rds("/volumes/USR1/yyan/project/tnbc_pre_atlas/export_to_paper/patient_meta_info/export.sample_info.rds")
    colnames(lib_patient_info)
    wanted_patient_info_list <- c("pCR_status", "archetype")

    lib_labid_to_response <- lib_patient_info[, c('lab_id', 'pCR_status')] %>% deframe()
    lib_ladid_to_archetype <- lib_patient_info[, c('lab_id', 'archetype')] %>% deframe()

}

#------ Create data without caring the patient groups info ------
cli_h2("Demultiplex into different object_suit regardless of any patient groups (pCR, archetypes)")

# object_suit <- object_suit_opts[14]
for (object_suit in object_suit_opts) {
    cli_h3(object_suit)
    dir_res <- file.path(dir_proj, object_suit, "inputs")
    fs::dir_create(dir_res)
    cell_idx_use <- NULL
    feature_idx_use <- NULL
    #------ demultiplex ------
    if (object_suit == "TME_only") {
        cell_idx_use <- df_meta$celltype %in% c("Mye", "T", "B", "Fibro", "Endo", "Peri")
        feature_idx_use <- colnames(df_value) %in% c(tme_feature_names, "cellname")
    } else if (object_suit == "Tumor_only") {
        cell_idx_use <- df_meta$celltype %in% c("Tumor")
        feature_idx_use <- colnames(df_value) %in% c(cancer_feature_names, "cellname")
    } else if (object_suit == "ALL") {
        cell_idx_use <- df_meta$celltype %in% c("Tumor", "Mye", "T", "B", "Fibro", "Endo", "Peri")
        feature_idx_use <- colnames(df_value) %in% c(cancer_feature_names, tme_feature_names, "cellname")
    } else {
        ## TME_and_*
        if (object_suit == "TME_and_tumor") {
            cell_idx_use <- df_meta$celltype %in% c("Tumor", "Mye", "T", "B", "Fibro", "Endo", "Peri")
            feature_idx_use <- colnames(df_value) %in% c(cancer_feature_names, tme_feature_names, "cellname")
        } else {
            ## TME_and_[*MP*]
            cancer_mp <- str_remove(object_suit, "TME_and_")
            if (!cancer_mp %in% cancer_feature_names) {
                warning(c("Invalid object_suit: ", object_suit))
                next()
            }
            cell_idx_use <- df_meta$celltype %in% c("Tumor", "Mye", "T", "B", "Fibro", "Endo", "Peri")
            feature_idx_use <- colnames(df_value) %in% c(cancer_mp, tme_feature_names, "cellname")
        }
    }

    print(table(cell_idx_use))
    df_meta_use <- df_meta[cell_idx_use, ]
    df_value_use <- df_value[cell_idx_use, ]
    df_value_use <- df_value_use[, feature_idx_use]

    if (object_suit == "TME_and_tumor") {
        ## merge the MPs into one column 'Tumor'
        df_value_use$Tumor <- as.numeric(rowSums(df_value_use[, cancer_feature_names, drop = FALSE], na.rm = TRUE) > 0)
        df_value_use <- df_value_use[, c("cellname", "Tumor", tme_feature_names)]
    }

    bad_features <- colnames(df_value_use)[sapply(df_value_use, function(x) sum(is.na(x))) != 0]
    if (length(bad_features) > 0) {
        cli_alert_warning("Check these features which have NA values in some cells")
        cli_ol(bad_features)
    }
    for (x in bad_features) {
        cli_alert_info(glue("Replacing NA values in {x} with 0"))
        df_value_use[[x]] <- replace_na(df_value_use[[x]], 0)
        rm(x)
    }

    stopifnot(identical(rownames(df_meta_use), rownames(df_value_use)))
    write_rds(df_meta_use, file.path(dir_res, "dataframe.cellmeta.rds"))
    write_rds(df_value_use, file.path(dir_res, "dataframe.value.rds"))

    ## write out the combo data frame
    feature_lvs <- colnames(select(df_value_use, -cellname))
    cell_assign <- select(df_value_use, -cellname) %>% apply(., 1, function(v) {
        idx <- which(v == 1)
        if (length(idx) == 1) {
            return(feature_lvs[idx])
        } else if (length(idx) == 0) {
            return("Unknown")
        } else {
            idx <- sample(idx, 1)
            return(feature_lvs[idx])
        }
    })
    df_combo_use <- df_meta_use
    stopifnot(identical(rownames(df_combo_use), names(cell_assign)))
    df_combo_use[[ident_is_what]] <- cell_assign
    df_combo_use <- df_combo_use[df_combo_use[[ident_is_what]] != "Unknown", ] # remove Unknown
    df_combo_use$pCR_status <- lib_labid_to_response[df_combo_use$sample]
    df_combo_use$archetype <- lib_ladid_to_archetype[df_combo_use$sample]
    write_rds(df_combo_use, file.path(dir_res, "dataframe.cellmeta_combo.rds"))

    #------ visualization ------
    set.seed(1026)
    df_meta_use_shuffle_ri <- sample(1:nrow(df_meta_use), nrow(df_meta_use)) ## shuffle cells
    #------ original spatial coords ------
    p1 <- ggplot(df_meta_use[df_meta_use_shuffle_ri, ], aes(x = coord_x, y = coord_y, color = sample)) +
        # geom_point_rast(size = .1) +
        geom_scattermore() +
        coord_equal() +
        scale_color_manual(values = pal_sample)
    ggsave(file.path(dir_res, "imgdimplot.adj_coord.1before.pdf"), p1, width = 10, height = 10, useDingbats = F)

    #------ first, adjust each tissue to the origin ------
    p2 <- ggplot(df_meta_use[df_meta_use_shuffle_ri, ], aes(x = coord_x_adj_to_origin, y = coord_y_adj_to_origin, color = sample)) +
        geom_scattermore() +
        # geom_point_rast(size = .1) +
        coord_equal() +
        scale_color_manual(values = pal_sample)
    ggsave(file.path(dir_res, "imgdimplot.adj_coord.2to_origin.pdf"), p2, width = 10, height = 10, useDingbats = F)

    #------ then, horitonally spread the samples ------
    p3 <- ggplot(df_meta_use[df_meta_use_shuffle_ri, ], aes(x = coord_x_adj, y = coord_y_adj, color = sample)) +
        # geom_point_rast(size = .1) +
        geom_scattermore() +
        coord_equal() +
        scale_color_manual(values = pal_sample)
    ggsave(file.path(dir_res, "imgdimplot.adj_coord.3after.pdf"),
        p3 + rremove("legend"),
        width = 3 * length(pal_sample), height = 4, useDingbats = F, limitsize = F
    )
    #------ Viz locations of cancer MP ------
    if (F) {
        # slow
        cli_alert_info("Visualizing locations of cancer MPs")
        df_cancer_mp <- df_combo_use %>% dplyr::filter(celltype %in% "Tumor")
        id_lvs <- unique(df_cancer_mp[[ident_is_what]])
        samples_to_collect <- unique(df_cancer_mp$sample)
        fs::dir_create(file.path(dir_res, "imgdimplot.cancer_idents"))
        for (sample_each in samples_to_collect) {
            cat("\n[", sample_each, "] ")
            df_sample_each <- df_combo_use %>% dplyr::filter(sample == sample_each)
            img_pdf_width_inch <- diff(range(df_sample_each$coord_x)) / 1000
            img_pdf_height_inch <- diff(range(df_sample_each$coord_y)) / 1000
            # img_pdf_width_inch <- max(img_pdf_width_inch, 5)
            # img_pdf_height_inch <- max(img_pdf_height_inch, 5)
            # img_pt_size <- max(diff(range(df_sample_each$coord_x)), diff(range(df_sample_each$coord_y))) / 1000/100 # no usage
            img_pt_size <- .5
            img_pt_scale <- min(diff(range(df_sample_each$coord_x)), diff(range(df_sample_each$coord_y))) / 1000 / 10
            for (id_lvs_each in id_lvs) {
                cat(sprintf(" %s", id_lvs_each))
                df_xx <- df_combo_use %>% dplyr::filter(
                    sample == sample_each,
                    celltype %in% "Tumor",
                    !!sym(ident_is_what) %in% id_lvs_each
                )
                pi <- ggplot(df_sample_each, aes(x = coord_x, y = coord_y)) +
                    # geom_scattermore() +
                    geom_point_rast(size = img_pt_size, color = "lightgrey", scale = img_pt_scale) +
                    geom_point_rast(
                        data = df_xx, aes(x = coord_x, y = coord_y),
                        size = img_pt_size, color = "red", scale = img_pt_scale
                    ) +
                    coord_equal() +
                    labs(title = sprintf("%s %s", sample_each, id_lvs_each)) +
                    theme(legend.position = "none") +
                    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
                    rremove("x.axis") +
                    rremove("y.axis")
                ggsave(
                    filename = file.path(
                        dir_res, "imgdimplot.cancer_idents",
                        sprintf("imgdimplot.%s_%s.pdf", sample_each, id_lvs_each)
                    ),
                    plot = pi,
                    width = img_pdf_width_inch,
                    height = img_pdf_height_inch + 0.1, useDingbats = FALSE, limitsize = TRUE
                )
            }
        }
        cat("\n")
    }
}



cat("[done]\n")
timestamp()
