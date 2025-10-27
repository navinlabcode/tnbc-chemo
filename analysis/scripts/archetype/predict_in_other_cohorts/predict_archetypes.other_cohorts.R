suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Predict archetypes for any other cohorts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-08-05 2025-03-17
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
    library(RcppML)
})

cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
    study_name <- cmdargs[1]
    f_type <- cmdargs[2]
} else {
    ## Examples of the available dataset ## 
    study_name <- "artemis_pretx"
    f_type <- "psbulk_allcells"
    
    study_name <- "2024Shiao"
    f_type <- "TNBC_baseline_cancercells"
    
    study_name <- "METABRIC"
    f_type <- "chemoYes"
    
    study_name <- "ISPY990"
    f_type <- "arm_Ctr"
    
    study_name <- "SCANB"
    f_type <- "dUTP"

    study_name <- "BrighTNess"
    f_type <- "arm_PaclitaxelCarboplatin"
    f_type <- "arm_ALL"
    f_type <- "arm_CHEMO"
    f_type <- "arm_Paclitaxel"

    study_name <- "NCPark2020"
    f_type <- "TNBCPreTXChemo"
    f_type <- "TNBCPreTX"
}
#------------------ ~~~ Load lib ~~~ --------------------
#------ lib ------
lib_mat <- read_rds(file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas/",
    "rds_rna-integrate/pat102/",
    "lv01.aneuploidy_tri_type.aneuploid.pure5/",
    "psbulk/fastnmf",
    "input.matrix.rds"
))
lib_arcs <- read_rds(
    file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/",
        "rds_rna-integrate/pat102/",
        "lv01.aneuploidy_tri_type.aneuploid.pure5/",
        "psbulk/fastnmf", "rank4",
        "deliver.patient_best_nmf.rds"
    )
)
lib_arcs <- gsub(x = lib_arcs, pattern = "fNMF", replacement = "ARC")
arc_lvs <- c("ARC1", "ARC2", "ARC3", "ARC4")
stopifnot(identical(names(lib_arcs), colnames(lib_mat)))
#------------------ ~~~ Input ~~~ --------------------
f_mat <- file.path(
    "/volumes/USR1/yyan/shared/data",
    study_name, f_type,
    "mat.rds"
)
f_df <- file.path(
    "/volumes/USR1/yyan/shared/data",
    study_name, f_type,
    "clinical.df.rds"
)
response_str <- "response"
response_lvs <- c("R", "NR")

pat_str <- NULL
if (study_name %in% c("ISPY990")) {
    response_str <- "Response"
    response_lvs <- c("pCR", "RD")
    pat_str <- "Patient Identifier"
}
if (study_name %in% c("BrighTNess")) {
    response_str <- "Response"
    response_lvs <- c("pCR", "RD")
    pat_str <- "patient_id"
}
if (study_name == "NCPark2020") {
    response_str <- "Response"
    response_lvs <- c("pCR", "RD")
    pat_str <- "sample_id"
}
if (study_name == "METABRIC") {
    # metabric's rest dataset
    f_mat <- file.path("/volumes/USR1/yyan/shared/data/METABRIC/ruli/", "tnbc_chemoYES_expression.mat.rds")
    f_df <- file.path("/volumes/USR1/yyan/shared/data/METABRIC/ruli/", "tnbc_chemoYES_survival.df.rds")
    pat_str <- "METABRIC_ID"
}
if (study_name == "SCANB") {
    # SCANB TNBC updated 2025-03-04
    if (f_type != "ALL_protocols") {
        dir_lib <- file.path("/volumes/USR1/yyan/shared/data/SCANB_TNBC_chemo/stringtie_fpkm_gene_data_unadjusted/", f_type)
    } else {
        dir_lib <- file.path("/volumes/USR1/yyan/shared/data/SCANB_TNBC_chemo/stringtie_fpkm_gene_data_LibProtocolAdjusted/", f_type)
    }
    f_df <- file.path(dir_lib, "scanb_tnbc.df_clinical.rds")
    f_mat <- file.path(dir_lib, "scanb_tnbc.expr_mat.rds")
    pat_str <- "GEX.assay"
}


if (study_name == "2021AyseBassez") {
    study_name <- "2021AyseBassez"
    f_type <- "TNBC_preTx_cancercells"
    f_mat <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("psbulk_normalized_DESeq2VST_matrix.%s.rds", f_type)
    )
    f_df <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("df_samplemeta_psbulk.%s.rds", f_type)
    )
    dir_res <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "archetypes_pred", f_type
    )
}

if (study_name == "2024Shiao") {
    study_name <- "2024Shiao"
    f_type <- "TNBC_baseline_cancercells"
    f_mat <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("psbulk_normalized_matrix.%s.rds", f_type)
    )
    f_df <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "data",
        sprintf("df_samplemeta_psbulk.%s.rds", f_type)
    )
    dir_res <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/other_public_cohorts",
        study_name, "archetypes_pred", f_type
    )
}

dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas/summary_fig",
    "archetype_prediction",
    study_name, f_type
)
fs::dir_create(dir_res)

cli_alert_info("Predicting archetypes for:\t {study_name} {f_type}")
cli_alert_info("Results will be saved to:\t {dir_res}")
#------------------ ~~~ Load input ~~~ --------------------
cli_h1("Load input")
mat_query <- read_rds(f_mat)
df_query <- read_rds(f_df)
if (!is.null(pat_str)) {
    df_query <- as.data.frame(df_query)
    rownames(df_query) <- df_query[[pat_str]]
}
df_query <- df_query[colnames(mat_query), ]
stopifnot(identical(colnames(mat_query), rownames(df_query)))



#------ center and remove negatives ------
mat_query <- t(mat_query)
mat_query <- scale(x = mat_query, center = TRUE, scale = FALSE)
mat_query <- t(mat_query)
mat_query[mat_query < 0] <- 0

shared_genes <- intersect(rownames(lib_mat), rownames(mat_query))
str(shared_genes)

#------------------ ~~~ Run NMF ~~~ --------------------
library(ComplexHeatmap)
mat <- cbind(lib_mat[shared_genes, ], mat_query[shared_genes, ])

if (any(is.na(mat))) {mat[is.na(mat)] <- 0}

res <- RcppML::nmf(mat, k = 4, tol = 1e-5, L1 = c(0.05, 0.05))
W <- res$w # genes x factor
H <- res$h # factor x samples
rownames(W) <- rownames(mat)
colnames(H) <- colnames(mat)
colnames(W) <- rownames(H) <- paste0("F", 1:4)
sample2nmf <- rownames(H)[apply(H, 2, nnet::which.is.max)]
names(sample2nmf) <- colnames(H)
head(sample2nmf)
write_rds(res, file.path(dir_res, "nmf_obj.rds"))

Heatmap(H[, colnames(lib_mat)], cluster_rows = F) %v%
    columnAnnotation(pred = sample2nmf[colnames(lib_mat)])
Heatmap(H[, colnames(mat_query)], cluster_rows = F) %v%
    columnAnnotation(pred = sample2nmf[colnames(mat_query)])

#------ Map the runtime factors to the archetypes ------
print(table(sample2nmf[names(lib_arcs)], lib_arcs))
write.csv(table(sample2nmf[names(lib_arcs)], lib_arcs), file.path(dir_res, "table.factor2ARCs.csv"))
absorb_labels <- function(query, ref) {
    tab <- as.data.frame.matrix(table(query, ref))
    ptab <- prop.table(tab)
    ref_l <- colnames(ptab)
    query_l <- rownames(ptab)
    j <- apply(ptab, 1, nnet::which.is.max)
    res <- structure(names = query_l, ref_l[j])
}
dict_factor2ARC <- absorb_labels(sample2nmf[names(lib_arcs)], lib_arcs)
dict_factor2ARC

#------ Decide archetypes for the input samples ------
arc_pred <- sample2nmf[colnames(mat_query)]
arc_pred <- structure(dict_factor2ARC[arc_pred],
    names = names(arc_pred)
)
table(arc_pred)
arc_pred <- factor(arc_pred, levels = arc_lvs)
write_csv(
    enframe(arc_pred, "sample_id", "archetype"),
    file.path(dir_res, "predicted_ARC.csv")
)
write_rds(arc_pred, file.path(dir_res, "predicted_ARC.rds"))

Heatmap(H[, colnames(mat_query)], cluster_rows = F) %v%
    columnAnnotation(pred = arc_pred[colnames(mat_query)])


#------------------ ~~~ Correctness on our training data ~~~ --------------------
cli_h1("Correctness on our training data")
arc_adhoc <- sample2nmf[colnames(lib_mat)]
arc_adhoc <- structure(dict_factor2ARC[arc_adhoc],
    names = names(arc_adhoc)
)
stopifnot(identical(names(arc_adhoc), names(lib_arcs)))
matched_rate <- sum(arc_adhoc == lib_arcs) / length(arc_adhoc)
write_lines(matched_rate, file.path(dir_res, "matched_rate.txt"))

message("Matched rate: ", signif(matched_rate, 3))

cat("Prediction of archetypes on any other cohort is finished. [done]\n\n")



cat("[done]")
timestamp()
