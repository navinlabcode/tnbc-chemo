#!/usr/bin/Rscript

##########################################
# Goal: Identify tumor cells using copykat
#
# Mix mode: Manually add normal breast tissue cells
#
# copykat_pred_default ## Ruli's default hclust=2
# copykat_pred_tirosh  ## CNV score v.s. CNV corr
# copykat_pred_leiden  ## clustering on iCNAs
#
##########################################
# Author: Yun Yan (yun.yan@uth.tmc.edu)
##########################################
library(docopt)
library(cli)
timestamp()
"Copykat_Mix Run copykat with an external reference tissue cells

Usage:
  copykat_mix.R <seurat3> <dir_copykat> [--name_sample=<abc>] [--downsample_demo] [--is_cellline]
  copykat_mix.R (-h | --help)
  copykat_mix.R --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  --name_sample=<abc>       Name for this sample [default: navinlab].
  --downsample_demo         If downsample to do demo.
  --is_cellline             If this is a cell line data.
" -> doc
arguments <- try(docopt(doc, version = "1.0.3"))
cli_rule("Running Arguments")
print(arguments)

suppressPackageStartupMessages({
    library(devtools)
    # load_all('/volumes//USR1/yyan/shared/pipelines/copykat')
    library(Seurat)
    library(copykat)
    library(tidyverse)
    library(reshape2)
    library(dplyr)
    library(purrr)
    library(ggplot2)
    library(ggpubr)
    library(RColorBrewer)
    theme_set(theme_pubr(base_size = 18, legend = "right"))
    library(ruok)
    library(glue)
    library(cli)
    library(tictoc)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(matrixStats)
    source("/volumes//USR1/yyan/project/tumor_plasticity/core/GTA/package_general_cna_analysis.R")
    source("/volumes//USR1/yyan/project/tumor_plasticity/core/pkg_add_rna_module_score.R")
    library(ComplexHeatmap)
    # ht_opt$fast_hclust <-  TRUE
})
fpath_input <- NULL

# TEST --------
# ## Example input
if (F) {
    fpath_input <- "/volumes/USR1/yyan/project/tnbc_visium_hd/data0/ART122/umi_count.matrix.rds"
    sample_name <- "ART122"
    dir_copykat <- file.path(
        "/volumes/USR1/yyan/project/tnbc_visium_hd/data",
        sample_name, "copykat_mix"
    )
}

## Specify copycat parameters [So far do not change]
feature_type <- "Symbol" # gene symbol / gene id
copykat_n_gene_per_chr_min <- 1
copykat_n_gene_per_seg_min <- 25
n_parallel <- 20
copykat_KS_cut <- 0.2
species_genome <- "hg38"
copykat_is_cellline <- "no"

# Run --------
if (try(any(class(arguments) != "try-error"))) {
    ## Inputs
    sample_name <- arguments$name_sample
    fpath_input <- arguments$seurat3
    ## Output
    dir_copykat <- arguments$dir_copykat
    fs::dir_create(dir_copykat)
    copykat_is_cellline <- ifelse(arguments$is_cellline, "yes", "no")
}

print(fpath_input)
print(sample_name)
print(c("is cell line: ", copykat_is_cellline))
print(dir_copykat)

## External reference cells
lib_infercnv_ref_sr3 <- read_rds(file.path(
    "/volumes//USR1/yyan/project/tumor_plasticity",
    "lib", "infercnv_ref", species_genome, "hbca_ref.sr3.rds"
)) ## Change if not breast cancer
# lib_infercnv_ref_sr3 <- UpdateSeuratObject(lib_infercnv_ref_sr3)
lib_gene_position <- readRDS(file.path(
    "/volumes//USR1/yyan/project/tumor_plasticity",
    "lib", "geneloc", glue("{species_genome}.ensembl_id.sorted.gr.rds")
))

fs::dir_create(dir_copykat)

## Let's rock
cli_rule("# Read Seurat object")
tic()
# sr3 <- readRDS(fpath_sr3);
# if this is a h5 file: 
if (grepl("\\.h5$", fpath_input)) {
    mat_obs <- Seurat::Read10X_h5(fpath_input)
} else if (grepl("\\.rds$", fpath_input)) {
    mat_obs <- readRDS(fpath_input)
} else {
    stop("Unknown input file type")
}
toc()

print(class(mat_obs))
if ("Seurat" %in% class(mat_obs)) {
    mat_obs <- GetAssayData(mat_obs, slot = "counts") %>% as.matrix()
}

# print(sr3)
# DefaultAssay(sr3) <- 'RNA'
print(dim(mat_obs)) ## genes x cells
obs_cellnames <- colnames(mat_obs)
str(obs_cellnames)
##
cli_rule("# Run copykat")
wd_orig <- getwd()
setwd(dir_copykat) ## copykat generates outputs silently
# mat_obs <- GetAssayData(object = sr3, slot = "counts", assay = 'RNA') %>% as.matrix()

mat_ref <- GetAssayData(lib_infercnv_ref_sr3, slot = "counts", assay = "RNA") %>% as.matrix()
cli_alert_info(sprintf("%d obs and %d ref cells are used to run copykat", ncol(mat_obs), ncol(mat_ref)))

idx <- match(rownames(mat_obs), rownames(mat_ref))
mat_ref <- mat_ref[idx, , drop = F]
colnames(mat_ref) <- paste0("ref", 1:ncol(lib_infercnv_ref_sr3))
mat_ref[is.na(mat_ref)] <- 0
rownames(mat_ref) <- rownames(mat_obs)

mat <- cbind(mat_obs, mat_ref)
# rm(mat_obs); rm(mat_ref)
class(mat)
mat <- as.matrix(mat)
## Remove chrX, chrY, chrM and non-standard chr genes
str(rownames(mat))
# lib_gene_position <- keepStandardChromosomes(lib_gene_position, pruning.mode = 'coarse')
# idx <- c(!as.character(seqnames(lib_gene_position)) %in% c('X', 'Y', 'MT'))
# idx <- rownames(mat) %in% lib_gene_position$symbol_seurat[idx]
# cat('Remove', sum(!idx), 'genes of input UMI are on chrX/Y/MT/non_std_chrom.\n')
# mat <- mat[which(idx), , drop=F]

fpath_tosave <- file.path(dir_copykat, "copykat_mix.rds")
if (!file.exists(fpath_tosave)) {
    tic()
    cpkt <- copykat(
        mat,
        id.type = feature_type,
        cell.line = copykat_is_cellline,
        ngene.chr = copykat_n_gene_per_chr_min,
        win.size = copykat_n_gene_per_seg_min,
        KS.cut = copykat_KS_cut,
        sam.name = sample_name,
        plot.genes = FALSE,
        n.cores = n_parallel
    )
    toc()
    cli_rule("# Export copykat result")
    saveRDS(object = cpkt, fpath_tosave)
} else {
    cli_alert_success(paste0("Load ", fpath_tosave))
    cpkt <- readRDS(fpath_tosave)
}
setwd(wd_orig)
# cpkt <- readRDS(fpath_tosave)
cli_alert_success(fpath_tosave)

combo_res <- data.frame(row.names = obs_cellnames, dummy=rep(NA, length(obs_cellnames)))
combo_res$Barcode <- obs_cellnames
##
cli_rule("# Visualize copykat prediction")
cpkt_tumor_pred <- as.data.frame(cpkt$prediction)
head(cpkt_tumor_pred)

## copykat seems have a bug in row names
if (T) {
    print(table(duplicated(cpkt_tumor_pred$cell.names)))
    bad_rows <- str_detect(rownames(cpkt_tumor_pred), "^X")
    print(table(bad_rows))
    cpkt_tumor_pred <- cpkt_tumor_pred[!bad_rows, ]
}
copykat_pred_default <- deframe(cpkt_tumor_pred)
copykat_pred_default <- copykat_pred_default[obs_cellnames]
copykat_pred_default <- replace_na(copykat_pred_default, 'Unknown')
combo_res$copykat_pred_default <- copykat_pred_default
print(table(combo_res$copykat_pred_default, useNA = 'ifany'))
rownames(combo_res) <- make.names(rownames(combo_res))

## export to Loupe browser
df_to_loupe <- cpkt_tumor_pred
colnames(df_to_loupe) <- c("Barcode", "copykat.pred")
write_csv(
    unique(df_to_loupe[df_to_loupe$Barcode %in% obs_cellnames, ]),
    file.path(dir_copykat, sprintf("%s.for_loupe.copykat_prediction.csv", sample_name))
)

if (F) { ## to take out
    idx <- match(colnames(sr3), cpkt_tumor_pred$cell.names)
    copykat_pred_default <- as.character(cpkt_tumor_pred$copykat.pred[idx])
    copykat_pred_default[is.na(copykat_pred_default)] <- "n/a"
    copykat_pred_default <- as.factor(copykat_pred_default)
    table(copykat_pred_default, useNA = "always")
    names(copykat_pred_default) <- Cells(sr3)
    sr3[["copykat_pred_default"]] <- copykat_pred_default ### --> to be exported
    table(sr3[["copykat_pred_default"]])
    color_aneuploid <- RColorBrewer::brewer.pal("Set1", n = 4)[4]
    color_diploid <- RColorBrewer::brewer.pal("Set1", n = 4)[3]
    color_copykat_ploidy <- c(
        `aneuploid` = color_aneuploid,
        `diploid` = color_diploid,
        `n/a` = "lightgrey"
    )

    p_scatter_dr <- DimPlot(
        sr3,
        reduction = "umap", group.by = "copykat_pred_default", order = "aneuploid"
    ) +
        theme_void() +
        scale_color_manual(
            values = color_copykat_ploidy,
            labels = pretty_table2str(table(sr3$copykat_pred_default))
        ) +
        coord_fixed() + labs(title = sample_name) + theme(legend.position = "bottom")

    ggsave(file.path(dir_copykat, "copykat_prediction-umap.png"),
        plot = p_scatter_dr, width = 4.2, height = 4
    )
    ggsave(file.path(dir_copykat, "copykat_prediction-umap.pdf"),
        plot = p_scatter_dr, width = 4.2, height = 4
    )
}

#-------------------------- Tirosh-based Prediction --------------------------
cli_rule("Tirosh Prediction")
decide_aneuploid <- function(X, idx_ref, idx_obs) {
    # X: matrix of bins x cells
    # idx_ref: index of the reference cells
    # idx_obs: index of the query cells
    # Yun Yan (yun.yan@uth.tmc.edu)
    mat_obs <- X[, idx_obs, drop = F]
    mat_ref <- X[, idx_ref, drop = F]
    #--------------------------
    # CNV score
    #--------------------------
    cnv_score_ref <- sqrt(colSums(mat_ref^2, na.rm = T))
    cnv_score_obs <- sqrt(colSums(mat_obs^2, na.rm = T))

    cnv_score_norm_window <- quantile(cnv_score_ref, c(0.5 / 100, 99.5 / 100))
    cnv_score_ref <- (cnv_score_ref - cnv_score_norm_window[1]) / diff(cnv_score_norm_window)
    cnv_score_obs <- (cnv_score_obs - cnv_score_norm_window[1]) / diff(cnv_score_norm_window)

    #--------------------------
    # CNV correlation
    #--------------------------
    mat_ref[mat_ref == log2(1)] <- log2(1 + 1e-4)
    mat_obs[mat_obs == log2(1)] <- log2(1 + 1e-4)
    # top_frac_model <- 10/100
    top_frac_model <- 1 / 100 ## Change June 1 2020

    idx_top_ref <- order(cnv_score_ref, decreasing = F) %>%
        head(., ceiling(top_frac_model * ncol(mat_ref)))
    idx_top_obs <- order(cnv_score_obs, decreasing = T) %>%
        head(., ceiling(top_frac_model * ncol(mat_obs)))
    cnv_model_ref <- mat_ref[, idx_top_ref, drop = F] %>% rowMeans(., na.rm = T)
    cnv_model_obs <- mat_obs[, idx_top_obs, drop = F] %>% rowMeans(., na.rm = T)
    cnv_model_ref[is.na(cnv_model_ref)] <- 0
    cnv_model_obs[is.na(cnv_model_obs)] <- 0

    cor.test(cnv_model_ref, cnv_model_obs)

    cor_ref_intra <- apply(mat_ref, 2, function(v) {
        cor(v, cnv_model_ref, method = "p")
    })
    # range(cor_ref_intra)
    cor_ref_obs <- apply(mat_ref, 2, function(v) {
        cor(v, cnv_model_obs, method = "p")
    })
    cnvcor_on_ref <- apply(mat_obs, 2, function(v) {
        cor(v, cnv_model_ref, method = "p")
    })
    # range(cnvcor_on_ref)

    cnvcor_on_obs <- apply(mat_obs, 2, function(v) {
        cor(v, cnv_model_obs, method = "p")
    })
    # range(cnvcor_on_obs)

    #--------------------------
    # Decide aneuploid v.s. normal
    #--------------------------

    # cnv_score_threshold <- median(cnv_score_ref) + 3 * mad(cnv_score_ref)
    # cnvcor_threshold <- median(cor_ref_obs) + 3 * mad(cor_ref_obs)
    cnv_score_threshold <- quantile(cnv_score_ref, 99 / 100, na.rm = T)
    if (cnv_score_threshold < 1) {
        cnv_score_threshold <- 1
    } ## Changed 2020-06-21
    cnvcor_threshold <- quantile(cor_ref_obs, 99 / 100, na.rm = T)
    if (is.na(cnvcor_threshold) | is.null(cnvcor_threshold) | is_empty(cnvcor_threshold)) {
        cnvcor_threshold <- 0.5
    } ## Changed 2020-06-21

    is_obs_cancer <- cnv_score_obs > cnv_score_threshold & cnvcor_on_obs > cnvcor_threshold
    is_obs_cancer[is.na(is_obs_cancer)] <- FALSE
    table(is_obs_cancer)
    head(is_obs_cancer)
    return(list(
        `is_aneuploid` = is_obs_cancer,
        `cnv_score` = cnv_score_obs,
        `cnv_cor` = cnvcor_on_obs,
        `cnv_score_ref` = cnv_score_ref,
        `cnv_cor_ref` = cor_ref_obs,
        `cnv_score_threshold` = cnv_score_threshold,
        `cnv_cor_threshold` = cnvcor_threshold
    ))
}
mat_cpkt <- cpkt$CNAmat[, 4:ncol(cpkt$CNAmat)] ## logCNV: bins x cells
mat_cpkt <- as.matrix(mat_cpkt) # bins x cells
print(table(cpkt_tumor_pred$copykat.pred))

stopifnot(identical(make.names(rownames(cpkt_tumor_pred)), colnames(mat_cpkt)))
colnames(mat_cpkt) <- rownames(cpkt_tumor_pred)
idx_obs <- match(intersect(make.names(obs_cellnames), colnames(mat_cpkt)), colnames(mat_cpkt))
idx_ref <- grep(pattern = "^ref", x = colnames(mat_cpkt))


prediction_copykat_tirosh <- decide_aneuploid(
    mat_cpkt, idx_ref, idx_obs
)
# prediction_copykat_tirosh <- readRDS(file.path(dir_copykat, 'copykat_aneuploidy_prediction_tirosh_intermediates.rds'))
pdf(file.path(dir_copykat, "copykat_aneuploidy_prediction_tirosh.pdf"), 4, 4, useDingbats = F)
plot(prediction_copykat_tirosh$cnv_score,
    prediction_copykat_tirosh$cnv_cor,
    col = prediction_copykat_tirosh$is_aneuploid + 1,
    xlim = range(c(prediction_copykat_tirosh$cnv_score, prediction_copykat_tirosh$cnv_score_ref)),
    ylim = c(min(c(prediction_copykat_tirosh$cnv_cor, prediction_copykat_tirosh$cnv_cor_ref)), 1),
    xlab = "CNV Score", ylab = "CNV Correlation",
    main = sprintf("Predicted aneuploid cells: %d", sum(prediction_copykat_tirosh$is_aneuploid))
)
points(prediction_copykat_tirosh$cnv_score_ref,
    prediction_copykat_tirosh$cnv_cor_ref,
    col = "blue"
)
abline(h = prediction_copykat_tirosh$cnv_cor_threshold, lty = "dashed")
abline(v = prediction_copykat_tirosh$cnv_score_threshold, lty = "dashed")
legend("bottomright",
    col = c("blue", "black", "red"),
    legend = c("Ref", "Obs-diploid", "Obs-tumor"), pch = 1
)
dev.off()

saveRDS(
    prediction_copykat_tirosh,
    file.path(dir_copykat, "copykat_aneuploidy_prediction_tirosh_intermediates.rds")
)

copykat_pred_tirosh <- ifelse(prediction_copykat_tirosh$is_aneuploid, "aneuploid", "diploid")
idx <- match(make.names(obs_cellnames), names(copykat_pred_tirosh))
copykat_pred_tirosh <- as.character(copykat_pred_tirosh[idx])
names(copykat_pred_tirosh) <- obs_cellnames
copykat_pred_tirosh <- replace_na(copykat_pred_tirosh, "Unknown")
copykat_pred_tirosh <- as.factor(copykat_pred_tirosh)
table(copykat_pred_tirosh)


# sr3[["copykat_pred_tirosh"]] <- copykat_pred_tirosh ### --> to be exported
combo_res$copykat_pred_tirosh <- copykat_pred_tirosh
print(table(combo_res$copykat_pred_tirosh))

if (F) { # to take out
  p_scatter_dr <- DimPlot(
    sr3,
    reduction = "umap", group.by = "copykat_pred_tirosh", order = "aneuploid"
  ) +
    theme_void() +
    scale_color_manual(
      values = color_copykat_ploidy,
      labels = pretty_table2str(table(sr3$copykat_pred_tirosh))
    ) +
    coord_fixed() + labs(title = sample_name) + theme(legend.position = "bottom")
  
  ggsave(file.path(dir_copykat, "copykat_prediction_tirosh-umap.png"),
         plot = p_scatter_dr, width = 4.2, height = 4
  )
  ggsave(file.path(dir_copykat, "copykat_prediction_tirosh-umap.pdf"),
         plot = p_scatter_dr, width = 4.2, height = 4
  )
}  
#-------------------------- iCNAs-clustering Prediction  --------------------------
cli_rule("Leiden iCNA prediction")
gta_cellanno <- combo_res[, c("copykat_pred_tirosh", "copykat_pred_default")]
gta_cellanno$category <- "obs"
# gta_cellanno$seurat_cluster <- Idents(sr3)
gta_cellanno$seurat_cluster <- NA

gta_ref_cellanno <- data.frame(row.names = colnames(mat_cpkt)[idx_ref])
gta_ref_cellanno$copykat_pred_tirosh <- "diploid"
gta_ref_cellanno$copykat_pred_default <- "diploid"
gta_ref_cellanno$category <- "ref"
gta_ref_cellanno$seurat_cluster <- NA

gta_cellanno <- rbind(gta_cellanno, gta_ref_cellanno)
gta_cellanno <- gta_cellanno[match(colnames(mat_cpkt), rownames(gta_cellanno)), ]
rownames(gta_cellanno) <- colnames(mat_cpkt)

gta_rowrange <- GRanges(cpkt$CNAmat[, 1], IRanges(cpkt$CNAmat[, 2], cpkt$CNAmat[, 2] + 1))

gta_mix <- SingleCellExperiment(
    assays = list(logCNA = mat_cpkt, CNA = 2^mat_cpkt - 1e-3),
    colData = gta_cellanno
)
rowRanges(gta_mix) <- gta_rowrange
gta_mix <- addPCA(gta_mix, npcs = 30, slot = "logCNA") # logCNA for ARTC02, ARTC94
gta_mix <- addUMAP(gta_mix, use_pcs = 1:30)
gta_mix <- run_phenograph(gta_mix, k = 50, use_dims = 1:2, reduction = "UMAP")
leiden_res <- 0.1
while (TRUE) {
    gta_mix <- run_leiden(gta_mix, resolution = leiden_res)
    leiden_str <- sprintf("Cluster_leiden_res_%s", leiden_res)
    leiden_n_cluster <- length(levels(colData(gta_mix)[[leiden_str]]))
    colData(gta_mix)[["icna_clusters"]] <- colData(gta_mix)[[leiden_str]]
    if (leiden_n_cluster >= 5) {
        break()
    } ## Maybe 1 obs-normal, 1 ref-normal, 2 obs-tumors, 1 obs-tumor/normal
    if (leiden_res >= 2) {
        break()
    } ## Higher resolutin won't help anymore; probably no aneuploidy cells at all.
    colData(gta_mix)[[leiden_str]] <- NULL
    leiden_res <- leiden_res + 0.1
}


p <- drplot(gta_mix, "UMAP", group.by = leiden_str) +
    scale_color_discrete(
        labels = pretty_table2str(table(colData(gta_mix)[[leiden_str]]))
    )
ggsave2(file.path(dir_copykat, "copykat_DNA_DR.leiden_cluster"),
    p,
    width = 6, height = 5.5
)
p <- drplot(gta_mix, "UMAP", group.by = "category", add_label = F)
ggsave2(file.path(dir_copykat, "copykat_DNA_DR.category"),
    p,
    width = 6, height = 5.5
)
# p
table(
    colData(gta_mix)[[leiden_str]],
    colData(gta_mix)[["copykat_pred_tirosh"]], useNA='ifany'
)

cli_alert_info("Find clusters without any reference cell; this cluster is probably the cancer cells.")
print(table(
    colData(gta_mix)[[leiden_str]],
    colData(gta_mix)[["category"]]
))
tmp <- table(
    colData(gta_mix)[[leiden_str]],
    colData(gta_mix)[["category"]]
) %>% as.matrix()
tmp_prop <- table(
    colData(gta_mix)[[leiden_str]],
    colData(gta_mix)[["category"]]
) %>%
    prop.table(., margin = 1) %>%
    as.matrix()

print(tmp)
leiden_is_tumor <- tmp[, "ref"] < 10 & tmp_prop[, "ref"] <= 1 / 100
leiden_is_tumor <- ifelse(leiden_is_tumor, "aneuploid", "diploid")

copykat_pred_leiden <- colData(gta_mix)[[leiden_str]]
copykat_pred_leiden <- recode_factor(copykat_pred_leiden, !!!leiden_is_tumor)
table(copykat_pred_leiden, useNA='ifany')
# copykat_pred_leiden <- factor(as.character(copykat_pred_leiden),
#     levels = c("aneuploid", "diploid")
# )
colData(gta_mix)[["copykat_pred_leiden"]] <- copykat_pred_leiden
idx <- match(make.names(obs_cellnames), colnames(gta_mix))
tmp <- as.character(colData(gta_mix)[["copykat_pred_leiden"]][idx])
tmp <- replace_na(tmp, "Unknown")
tmp <- factor(as.character(tmp))#, levels = c("aneuploid", "diploid"))
# sr3[["copykat_pred_leiden"]] <- tmp ## --> to be exported
combo_res$copykat_pred_leiden <- tmp

table(
    colData(gta_mix)[[leiden_str]],
    colData(gta_mix)[["copykat_pred_leiden"]]
)

if (F) { # to take out
  p_scatter_dr <- DimPlot(
    sr3,
    reduction = "umap", group.by = "copykat_pred_leiden", order = "aneuploid"
  ) +
    theme_void() +
    scale_color_manual(
      values = color_copykat_ploidy,
      labels = pretty_table2str(table(sr3$copykat_pred_leiden))
    ) +
    coord_fixed() + labs(title = sample_name) + theme(legend.position = "bottom")
  ggsave(file.path(dir_copykat, "copykat_pred_leiden-umap.png"),
         plot = p_scatter_dr, width = 4.2, height = 4
  )
  ggsave(file.path(dir_copykat, "copykat_pred_leiden-umap.pdf"),
         plot = p_scatter_dr, width = 4.2, height = 4
  )
}

## Viz icna_leiden clusters
tmp <- colData(gta_mix[, intersect(make.names(obs_cellnames), colnames(gta_mix))])[["icna_clusters"]]
names(tmp) <- intersect(make.names(obs_cellnames), colnames(gta_mix))
tmp <- tmp[match(make.names(obs_cellnames), names(tmp))]
combo_res$icna_leiden <- tmp
if (F) { # to take out
  sr3$icna_leiden <- tmp
  p_scatter_dr <- UMAPPlot(sr3, group.by = "icna_leiden") + theme_void(base_size = 20) +
    labs(title = sample_name, color = "icna_leiden") +
    coord_equal() + theme(legend.position = "bottom")
  ggsave(file.path(dir_copykat, "copykat_icna_leiden-umap.png"),
         plot = p_scatter_dr, width = 4.2, height = 4
  )
  ggsave(file.path(dir_copykat, "copykat_icna_leiden-umap.pdf"),
         plot = p_scatter_dr, width = 4.2, height = 4
  )
}
rm(tmp)


#-------------------------- Export Results --------------------------
# copykat_pred_default ## Ruli's default hclust=2
# copykat_pred_tirosh  ## CNV score v.s. CNV corr
# copykat_pred_leiden  ## clustering on iCNAs
combo_res$dummy <- NULL
print(tail(combo_res))
write_tsv(combo_res, file.path(dir_copykat, "copykat_pred_report.tsv"))

combo_res_origcell <- combo_res
rownames(combo_res_origcell) <- combo_res_origcell$Barcode
stopifnot(identical(rownames(combo_res_origcell), obs_cellnames))
write_rds(combo_res_origcell, file.path(dir_copykat, "copykat_pred_report.seurat3_meta.rds"))

write_rds(gta_mix, file.path(dir_copykat, "copykat_mix.sce.rds"))

combo_res_origcell <- read_rds(file.path(dir_copykat, "copykat_pred_report.seurat3_meta.rds"))
for (z in c('copykat_pred_default', 'copykat_pred_tirosh', 
            'copykat_pred_leiden', 'icna_leiden')) {
  if (! z %in% colnames(combo_res_origcell)) {next()}
  combo_res_origcell[, c('Barcode', z)] %>%
    write_csv(., file.path(dir_copykat, sprintf('for_loupe.%s.%s.csv', sample_name, z)))
}
cli_rule("[DONE]")
timestamp()


#-------------------------- Final Visualization --------------------------
cli_rule("Final Heatmap Viz")
set.seed(22)

hm_obj_sc <- plot_heatmap_sc_manual(
  gta_mix,
  clip = c(-1, 1),
  cell_group_by = 'copykat_pred_tirosh',
  anno_rows_category = c(
    "category",
    "copykat_pred_default",
    "copykat_pred_tirosh",
    "copykat_pred_leiden"
  )
)

pdf(file.path(dir_copykat, "copykat_heatmap3.copykat_pred_tirosh.pdf"), width = 13, height = 10)
draw(hm_obj_sc,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom"
)
dev.off()

hm_obj_sc <- plot_heatmap_sc_manual(
  gta_mix,
  clip = c(-1, 1),
  cell_group_by = leiden_str,
  anno_rows_category = c(
    "category",
    "copykat_pred_default",
    "copykat_pred_tirosh",
    "copykat_pred_leiden"
  )
)
pdf(file.path(dir_copykat, "copykat_heatmap3.pdf"), width = 13, height = 10)
draw(hm_obj_sc,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom"
)
dev.off()

# DefaultAssay(sr3) <- 'RNA'
# feature_plot(sr3, G=c('FABP7', 'HPGD', 'ACSM3',
#                       'EPCAM', 'PTPRC', 'COL1A1'), dims = c('UMAP_1', 'UMAP_2'),
#              data_slot = 'scale.data',
#              cutoff.max = 3, cutoff.min = -3,
#              pt.alpha = 1, do.order = T,
#              ncol=6)
# feature_plot(sr3, G=c('EPCAM', 'PTPRC', 'COL1A1'), dims = c('UMAP_1', 'UMAP_2'),
#              data_slot = 'scale.data',
#              cutoff.max = 3, cutoff.min = -3,
#              pt.alpha = .5, do.order = F,
#              ncol=3)
# tmp <- colData(gta_mix[, intersect(Cells(sr3), colnames(gta_mix))])[[leiden_str]]
# names(tmp) <- intersect(Cells(sr3), colnames(gta_mix))
# tmp <- tmp[match(Cells(sr3), names(tmp))]
# sr3$icna_leiden <- tmp
# UMAPPlot(sr3, group.by='icna_leiden') + theme_void(base_size = 20)+
#   coord_equal() + theme(legend.position = 'bottom')

cat('[DONE] CopyKAT prediction and visualization\n')
timestamp()
