suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Run TransferLabel of Seurat v4 in the Evaluation mode
# Use 9 portions of genes to run TransferLabel and use the
# leftover 1 portion of genes to evaluate the performance of
# the determined Anchor.
#
#
# Known problem: If too many cells, the step 'TransferData - genes' fails. 
# Solution: split the query cells into smaller groups and then combinle the resulting matrix `imputation@assays$id@data`.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: text
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
    library(Seurat) # 4.3.0
    library(future)
    # options(future.globals.maxSize = 20 * 1024^3)
    options(future.globals.maxSize = 50 * 1024^3)
})
cmdargs <- commandArgs(trailingOnly = TRUE)
print(cmdargs)
param.FindTransferAnchors_n_dims <- 30
param.run_TransferData_genes <- FALSE
if (length(cmdargs) > 0) {
    f_sc <- cmdargs[[1]]
    f_sp <- cmdargs[[2]]
    cat_transfer <- cmdargs[[3]] # celltype | cell_state_paper
    sc_assay <- cmdargs[[4]] # RNA
    sp_assay <- cmdargs[[5]] # Xenium
    dir_gene_portions <- cmdargs[[6]] # xenium5k | xenium5kPlus
    gene_portion_index <- cmdargs[[7]] # 1,2,...10
    gene_portion_index <- as.numeric(gene_portion_index)
} else {
    if (T) {
        f_sc <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/downto_1000_by_cell_state_paper/ready.sr3.rds"
        cat_transfer <- "celltype"
        sc_assay <- "RNA"
    }
    f_sp <- "/volumes/USR1/yyan/project/tnbc_xenium/data/ART23/nonbinarized_pca/xenium_nonbinarized_pca.seurat.rds"
    sp_assay <- "Xenium"
    gene_portion_index <- 4
    dir_gene_portions <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/downto_1000_by_cell_state_paper/gene_portions_xenium5k"
}

dir_anchor <- file.path(
    dirname(f_sp),
    sprintf(
        "TransferLabelEval.ref_%s.query_%s.gene_portion_%s",
        sc_assay, sp_assay, gene_portion_index
    )
)
dir_cat_transfer <- file.path(dir_anchor, cat_transfer)
dir_anchor_eval <- file.path(dir_anchor, "eval_anchors")
fs::dir_create(dir_anchor)
fs::dir_create(dir_cat_transfer)
fs::dir_create(dir_anchor_eval)
#------------------ ~~~ Readin ~~~ --------------------
sc <- read_rds(f_sc)
sp <- read_rds(f_sp)

DefaultAssay(sc) <- sc_assay
DefaultAssay(sp) <- sp_assay

adhoc_load_gene_portions <- function(
    dir_gene_portions, eval_gene_portion_index, n_portions = 10) {
    f_test <- file.path(
        dir_gene_portions,
        sprintf("genes_of_portion_%d.rds", eval_gene_portion_index)
    )
    f_train_list <- file.path(
        dir_gene_portions,
        sprintf(
            "genes_of_portion_%d.rds",
            setdiff(1:n_portions, eval_gene_portion_index)
        )
    )
    g_test <- read_rds(f_test)
    g_train <- as.character(unlist(lapply(f_train_list, read_rds)))
    return(list(
        genes_test = g_test,
        genes_train = g_train
    ))
}

## read genes for training and testing
genes_combo <- adhoc_load_gene_portions(
    dir_gene_portions = dir_gene_portions,
    eval_gene_portion_index = gene_portion_index
)
str(genes_combo)

#------------------ ~~~ FindAnchor ~~~ --------------------

cli_h1("FindTransferAnchors")
timestamp()
f_anchors <- file.path(dir_anchor, "anchors.rds")
if (!file.exists(f_anchors)) {
    anchors <- FindTransferAnchors(
        reference = sc, query = sp,
        reference.assay = sc_assay,
        query.assay = sp_assay,
        features = genes_combo$genes_train,
        reduction = "cca",
        k.filter = NA,
        dims = 1:param.FindTransferAnchors_n_dims
    )
    write_rds(anchors, f_anchors)
    cli_alert_success("[finished] FindTransferAnchors.")
} else {
    cli_alert_success("reading the existing anchors...")
    anchors <- read_rds(f_anchors)
}
timestamp()

#------------------ ~~~ TransferData - genes ~~~ --------------------
# Ref: justifying 'weight.reduction': https://github.com/satijalab/seurat/issues/2926
## OPTIONAL: TransferData for genes
if (param.run_TransferData_genes) {
    cli_h1("TransferData - genes")
    timestamp()
    f_impute <- file.path(dir_anchor_eval, "obj.TransferData_impute_genes.rds")
    f_impute_data <- file.path(dir_anchor_eval, "TransferData_impute_genes.matrix.rds")
    if (!file.exists(f_impute)) {
        refdata <- GetAssayData(sc, slot = "data", assay = sc_assay)
        imputation <- TransferData(
            anchorset = anchors,
            query = sp,
            refdata = refdata,
            weight.reduction = "pca",
            dims = 1:param.FindTransferAnchors_n_dims,
            slot = "data"
        )
        write_rds(imputation, f_impute)
        cli_alert_success("[finished] TransferData")
    } else {
        cli_alert_success("reading the existing TransferData...")
        imputation <- read_rds(f_impute)
    }
    timestamp()
    class(imputation) # seurat
    print(imputation@assays$id) # the new/imputed assay gene data
    print(dim(imputation@assays$id@data))
    if (!file.exists(f_impute_data)) {
        impute_genes_all <- imputation@assays$id@data
        write_rds(imputation@assays$id@data, f_impute_data)
    } else {
        impute_genes_all <- read_rds(f_impute_data)
    }


    f_sp_impute_g_tain <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_imputed_genes_train"))
    f_sp_impute_g_test <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_imputed_genes_test"))
    f_sp_observe_g_train <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_observed_genes_train"))
    f_sp_observe_g_test <- file.path(dir_anchor_eval, sprintf("%s.matrix.rds", "sp_observed_genes_test"))

    if (!file_exists(f_sp_observe_g_test)) {
        ## impute_sp_g_train + obs_sp_g_train: evaluating training performance
        ## impute_sp_g_test + obs_sp_g_test  : evaluating testing performance

        impute_sp_g_train <- impute_genes_all[genes_combo$genes_train, ]
        impute_sp_g_test <- impute_genes_all[genes_combo$genes_test, ]

        obs_sp_g_train <- GetAssayData(sp, slot = "data", assay = sp_assay)[genes_combo$genes_train, ]
        obs_sp_g_test <- GetAssayData(sp, slot = "data", assay = sp_assay)[genes_combo$genes_test, ]
        class(impute_sp_g_train)
        class(impute_sp_g_test)
        class(obs_sp_g_train)
        class(obs_sp_g_test)

        stopifnot(identical(colnames(impute_sp_g_train), Cells(sp)))
        stopifnot(identical(colnames(impute_sp_g_test), Cells(sp)))
        stopifnot(identical(colnames(obs_sp_g_train), Cells(sp)))
        stopifnot(identical(colnames(obs_sp_g_test), Cells(sp)))

        write_rds(impute_sp_g_train, f_sp_impute_g_tain)
        write_rds(impute_sp_g_test, f_sp_impute_g_test)
        write_rds(obs_sp_g_train, f_sp_observe_g_train)
        write_rds(obs_sp_g_test, f_sp_observe_g_test)
    } else {
        impute_sp_g_train <- read_rds(f_sp_impute_g_tain)
        impute_sp_g_test <- read_rds(f_sp_impute_g_test)
        obs_sp_g_train <- read_rds(f_sp_observe_g_train)
        obs_sp_g_test <- read_rds(f_sp_observe_g_test)
    }
}
#------------------ ~~~ TransferData - categorical labels ~~~ --------------------
cli_h1(c("TransferData - ", cat_transfer))
timestamp()
f_labeltransfer <- file.path(dir_cat_transfer, sprintf("TransferData_%s.dataframe.rds", cat_transfer))
if (!file_exists(f_labeltransfer)) {
    pred <- TransferData(
        anchorset = anchors,
        query = sp,
        refdata = sc@meta.data[[cat_transfer]],
        weight.reduction = "pca",
        dims = 1:param.FindTransferAnchors_n_dims
    )
    pred_score_matrix <- t(pred@assays$prediction.score.id@data)
    identical(rownames(pred@meta.data), rownames(pred_score_matrix))

    colnames(pred_score_matrix) <- paste0("prediction.score.", colnames(pred_score_matrix))
    pred_id <- pred@meta.data$predicted.id
    pred_score_max <- pred@meta.data$predicted.id.score
    # all.equal.numeric(apply(pred_score_matrix, 1, max), pred_score_max)

    class(pred)
    colnames(pred@meta.data)
    table(pred$predicted.id)
    hist(pred$predicted.id.score)
    #
    # [1] "predicted.id"           "prediction.score.T"
    # [3] "prediction.score.Mye"   "prediction.score.B"
    # [5] "prediction.score.Fibro" "prediction.score.Endo"
    # [7] "prediction.score.Peri"  "prediction.score.Tumor"
    # [9] "prediction.score.max"

    pred_score_matrix <- as.data.frame(pred_score_matrix)
    pred_df <- cbind(
        predicted.id = pred$predicted.id,
        prediction.score.max = pred$predicted.id.score,
        pred_score_matrix
    )
    head(pred_df)

    write_rds(pred_df, f_labeltransfer)
    write_csv(pred_df, paste0(f_labeltransfer, ".csv"))


    if (FALSE) {
        pred <- TransferData(
            anchorset = anchors,
            refdata = sc@meta.data[[cat_transfer]],
            weight.reduction = "cca",
            dims = 1:param.FindTransferAnchors_n_dims
        )
        class(pred)
        colnames(pred)
        hist(pred$prediction.score.max)
        table(pred$predicted.id)
    }
} else {
    pred_df <- read_rds(f_labeltransfer)
}
timestamp()


#------------------ ~~~ Evaluate transfered genes ~~~ --------------------
# See 'xenium.celltype.TransferLabel_EvalMode.eval_gene.R'

#------------------ ~~~ Visualize - transferred genes ~~~ --------------------

#------------------ ~~~ Visualize - categorical labels ~~~ --------------------

#------------------ ~~~ Export transferred label to Xenium Explorer ~~~ --------------------
source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
for (ident_str in c("predicted.id")) {
    df_ident <- pred_df %>%
        rownames_to_column("cell_name") %>%
        df_to_exnium_explorer(., "cell_name", ident_str)
    write_csv(df_ident, file.path(
        dir_cat_transfer,
        sprintf("to_xenium_explorer_groups.%s.csv", ident_str)
    ))
}

cat("[done] R script")
timestamp()
