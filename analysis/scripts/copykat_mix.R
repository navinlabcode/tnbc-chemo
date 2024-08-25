#!/usr/bin/Rscript

##########################################
# Goal: Identify tumor cells using copykat
# 
# Mix mode: Manually add normal breast tissue cells
# 
# copykat_pred_default ## Ruli's default hclust=2
# copykat_pred_tirosh  ## CNV score v.s. CNV corr
# copykat_pred_leiden  ## clustering on iCNAs
# copykat_pred_livnat  ## corr with tumor/normal bulk RNA samples
# 
##########################################
# Author: Yun Yan (yun.yan@uth.tmc.edu)
##########################################
library(docopt); library(cli); timestamp()
'Copykat_Mix Run copykat with an external reference tissue cells

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
' -> doc
arguments <- try(docopt(doc, version = '1.0.3'))
cli_rule('Running Arguments')
print(arguments)

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(tidyverse); library(reshape2); library(dplyr); library(purrr)
  library(ggplot2); library(ggpubr); library(RColorBrewer)
  theme_set(theme_pubr(base_size=18, legend = 'right'))
  library(glue); library(cli); library(tictoc)
  library(GenomicRanges); library(GenomeInfoDb)
  library(matrixStats)
  library(ComplexHeatmap)
  # ht_opt$fast_hclust <-  TRUE
})

# TEST --------
# ## Inputs
dir_project <- 'data/applications/tnbc_chemo'
sample_name <- "art089_pre_core" 
fpath_sr3 <- glue('{dir_project}/rds_rna/{sample_name}/yansr3/ready.sr3.rds')
# # Output
dir_copykat <- glue('{dir_project}/rds_rna/{sample_name}/copykat_mix'); fs::dir_create(dir_copykat)


## Specify copycat parameters [So far do not change]
feature_type <- 'Symbol'  # gene symbol / gene id
copykat_n_gene_per_chr_min <- 1
copykat_n_gene_per_seg_min <- 25
n_parallel <- 20
copykat_KS_cut <- 0.2
species_genome <- 'hg38'
copykat_is_cellline <- 'no'

# Run --------
if (try(any(class(arguments) != 'try-error'))) {
  ## Inputs
  sample_name <- arguments$name_sample
  fpath_sr3 <- arguments$seurat3
  ## Output
  dir_copykat <- arguments$dir_copykat; fs::dir_create(dir_copykat)
  copykat_is_cellline <- ifelse(arguments$is_cellline, 'yes', 'no')
}

## External reference cells
lib_infercnv_ref_sr3 <- readRDS(file.path(
  "data", 
  "lib", "infercnv_ref", species_genome, 'hbca_ref.sr3.rds'))  ## Change if not breast cancer


## Let's rock
cli_rule('# Read Seurat object'); tic()
sr3 <- readRDS(fpath_sr3); toc()
print(sr3)
DefaultAssay(sr3) <- 'RNA'


##
cli_rule('# Run copykat')
wd_orig <- getwd(); setwd(dir_copykat) ## copykat generates outputs silently
mat_obs <- GetAssayData(object = sr3, slot = "counts", assay = 'RNA') %>% as.matrix()
mat_ref <- GetAssayData(lib_infercnv_ref_sr3, slot='counts', assay = 'RNA') %>% as.matrix()
cli_alert_info(sprintf('%d obs and %d ref cells are used to run copykat', ncol(mat_obs), ncol(mat_ref)))

idx <- match(rownames(mat_obs), rownames(mat_ref))
mat_ref <- mat_ref[idx, , drop=F]
colnames(mat_ref) <- paste0('ref', 1:ncol(lib_infercnv_ref_sr3))
mat_ref[is.na(mat_ref)] <- 0
rownames(mat_ref) <- rownames(mat_obs)

mat <- cbind(mat_obs, mat_ref)
rm(mat_obs); rm(mat_ref)


fpath_tosave <- file.path(dir_copykat, 'copykat_mix.rds')
if (! file.exists(fpath_tosave)) {
  tic()
  cpkt <- copykat(
    mat, id.type = feature_type, 
    cell.line = copykat_is_cellline,
    ngene.chr = copykat_n_gene_per_chr_min, 
    win.size = copykat_n_gene_per_seg_min, 
    KS.cut = copykat_KS_cut,
    sam.name = sample_name,
    n.cores = n_parallel)
  toc()
  cli_rule('# Export copykat result')
  saveRDS(object = cpkt, fpath_tosave)
} else {
  cli_alert_success(paste0('Load ', fpath_tosave))
  cpkt <- readRDS(fpath_tosave)
}
setwd(wd_orig)
# cpkt <- readRDS(fpath_tosave)
cli_alert_success(fpath_tosave)

##
cli_rule('# Visualize copykat prediction')
cpkt_tumor_pred <- as.data.frame(cpkt$prediction)
head(cpkt_tumor_pred)

idx <- match(colnames(sr3), cpkt_tumor_pred$cell.names)
copykat_pred_default <- as.character(cpkt_tumor_pred$copykat.pred[idx])
copykat_pred_default[is.na(copykat_pred_default)] <- 'n/a'
copykat_pred_default <- as.factor(copykat_pred_default)
table(copykat_pred_default, useNA='always')
names(copykat_pred_default) <- Cells(sr3)
sr3[['copykat_pred_default']] <- copykat_pred_default   ###--> to be exported
table(sr3[['copykat_pred_default']])
color_aneuploid <- RColorBrewer::brewer.pal('Set1', n=4)[4]
color_diploid <- RColorBrewer::brewer.pal('Set1', n=4)[3]
color_copykat_ploidy <- c(`aneuploid`=color_aneuploid, 
                          `diploid`=color_diploid,
                          `n/a`='lightgrey')

p_scatter_dr <- DimPlot(
  sr3, reduction = 'umap', group.by='copykat_pred_default', order = 'aneuploid') +
  theme_void() + 
  scale_color_manual(values = color_copykat_ploidy,
                     labels = pretty_table2str(table(sr3$copykat_pred_default))) +
  coord_fixed() + labs(title=sample_name) + theme(legend.position = 'bottom')

ggsave(file.path(dir_copykat, 'copykat_prediction-umap.pdf'), 
       plot = p_scatter_dr, width = 4.2, height = 4)     


#-------------------------- Tirosh-based Prediction --------------------------  
cli_rule('Tirosh Prediction')
decide_aneuploid <- function(X, idx_ref, idx_obs){
  # X: matrix of bins x cells
  # idx_ref: index of the reference cells 
  # idx_obs: index of the query cells
  # Yun Yan (yun.yan@uth.tmc.edu)
  mat_obs <- X[, idx_obs, drop=F]
  mat_ref <- X[, idx_ref, drop=F]
  #--------------------------
  # CNV score
  #--------------------------
  cnv_score_ref <- sqrt(colSums(mat_ref ^ 2, na.rm = T))
  cnv_score_obs <- sqrt(colSums(mat_obs ^ 2, na.rm = T))
  
  cnv_score_norm_window <- quantile(cnv_score_ref, c(0.5/100, 99.5/100))
  cnv_score_ref <- (cnv_score_ref - cnv_score_norm_window[1]) / diff(cnv_score_norm_window)
  cnv_score_obs <- (cnv_score_obs - cnv_score_norm_window[1]) / diff(cnv_score_norm_window)
  
  #--------------------------
  # CNV correlation
  #--------------------------
  mat_ref[mat_ref == log2(1)] <- log2(1+1e-4)
  mat_obs[mat_obs == log2(1)] <- log2(1+1e-4)

  top_frac_model <- 1/100 
  
  idx_top_ref <- order(cnv_score_ref, decreasing = F) %>%
    head(., ceiling(top_frac_model * ncol(mat_ref)))
  idx_top_obs <- order(cnv_score_obs, decreasing = T) %>%
    head(., ceiling(top_frac_model * ncol(mat_obs)))
  cnv_model_ref <- mat_ref[, idx_top_ref, drop=F] %>% rowMeans(., na.rm = T)
  cnv_model_obs <- mat_obs[, idx_top_obs, drop=F] %>% rowMeans(., na.rm = T)
  cnv_model_ref[is.na(cnv_model_ref)] <- 0
  cnv_model_obs[is.na(cnv_model_obs)] <- 0
  
  cor.test(cnv_model_ref, cnv_model_obs)
  
  cor_ref_intra <- apply(mat_ref, 2, function(v){
    cor(v, cnv_model_ref, method = 'p')
  })
  # range(cor_ref_intra)
  cor_ref_obs <- apply(mat_ref, 2, function(v){
    cor(v, cnv_model_obs, method='p')
  })
  cnvcor_on_ref <- apply(mat_obs, 2, function(v){
    cor(v, cnv_model_ref, method = 'p')
  })
  # range(cnvcor_on_ref)
  
  cnvcor_on_obs <- apply(mat_obs, 2, function(v){
    cor(v, cnv_model_obs, method='p')
  })
  # range(cnvcor_on_obs)
  
  #--------------------------
  # Decide aneuploid v.s. normal
  #--------------------------
  
  cnv_score_threshold <- quantile(cnv_score_ref, 99/100, na.rm = T)
  if (cnv_score_threshold < 1) {cnv_score_threshold <- 1}  
  cnvcor_threshold <- quantile(cor_ref_obs, 99/100, na.rm = T)
  if (is.na(cnvcor_threshold) | is.null(cnvcor_threshold) | is_empty(cnvcor_threshold)) {
    cnvcor_threshold <- 0.5 
  } 

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
mat_cpkt <- cpkt$CNAmat[, 4:ncol(cpkt$CNAmat)]  ## logCNV: bins x cells 
mat_cpkt <- as.matrix(mat_cpkt)
stopifnot(identical(make.names(rownames(cpkt_tumor_pred)), colnames(mat_cpkt)))
colnames(mat_cpkt) <- rownames(cpkt_tumor_pred)
idx_obs <- match(intersect(Cells(sr3), colnames(mat_cpkt)), colnames(mat_cpkt))
idx_ref <- grep(pattern = '^ref', x = colnames(mat_cpkt))


prediction_copykat_tirosh <- decide_aneuploid(
  mat_cpkt, idx_ref, idx_obs)

pdf(file.path(dir_copykat, 'copykat_aneuploidy_prediction_tirosh.pdf'), 4, 4, useDingbats = F)
plot(prediction_copykat_tirosh$cnv_score, 
     prediction_copykat_tirosh$cnv_cor,
     col=prediction_copykat_tirosh$is_aneuploid+1,
     xlim=range(c(prediction_copykat_tirosh$cnv_score, prediction_copykat_tirosh$cnv_score_ref)),
     ylim=c(min(c(prediction_copykat_tirosh$cnv_cor, prediction_copykat_tirosh$cnv_cor_ref)),1),
     xlab='CNV Score', ylab='CNV Correlation',
     main=sprintf('Predicted aneuploid cells: %d', sum(prediction_copykat_tirosh$is_aneuploid)))
points(prediction_copykat_tirosh$cnv_score_ref, 
       prediction_copykat_tirosh$cnv_cor_ref, col='blue')
abline(h=prediction_copykat_tirosh$cnv_cor_threshold, lty='dashed')
abline(v=prediction_copykat_tirosh$cnv_score_threshold, lty='dashed')
legend('bottomright', col=c('blue', 'black', 'red'), 
       legend=c('Ref', 'Obs-diploid', 'Obs-tumor'), pch = 1)
dev.off()

saveRDS(prediction_copykat_tirosh, 
        file.path(dir_copykat, 'copykat_aneuploidy_prediction_tirosh_intermediates.rds'))

copykat_pred_tirosh <- ifelse(prediction_copykat_tirosh$is_aneuploid, 'aneuploid', 'diploid')
idx <- match(Cells(sr3), names(copykat_pred_tirosh))
copykat_pred_tirosh <- as.character(copykat_pred_tirosh[idx])
names(copykat_pred_tirosh) <- Cells(sr3)
copykat_pred_tirosh <- replace_na(copykat_pred_tirosh, 'diploid')
copykat_pred_tirosh <- as.factor(copykat_pred_tirosh)
table(copykat_pred_tirosh)

sr3[['copykat_pred_tirosh']] <- copykat_pred_tirosh   ###--> to be exported
print(table(sr3$copykat_pred_tirosh))
p_scatter_dr <- DimPlot(
  sr3, reduction = 'umap', group.by='copykat_pred_tirosh', order = 'aneuploid') +
  theme_void() + 
  scale_color_manual(values = color_copykat_ploidy,
                     labels = pretty_table2str(table(sr3$copykat_pred_tirosh))) +
  coord_fixed() + labs(title=sample_name) + theme(legend.position = 'bottom')

ggsave(file.path(dir_copykat, 'copykat_prediction_tirosh-umap.pdf'), 
       plot = p_scatter_dr, width = 4.2, height = 4)  

#-------------------------- iCNAs-clustering Prediction  --------------------------  
cli_rule('Leiden iCNA prediction')
gta_cellanno <- sr3[[]][, c('copykat_pred_tirosh', 'copykat_pred_default')]
gta_cellanno$category <- 'obs'
gta_cellanno$seurat_cluster <- Idents(sr3)

gta_ref_cellanno <- data.frame(row.names = colnames(mat_cpkt)[idx_ref])
gta_ref_cellanno$copykat_pred_tirosh <- 'diploid'
gta_ref_cellanno$copykat_pred_default <- 'diploid'
gta_ref_cellanno$category <- 'ref'
gta_ref_cellanno$seurat_cluster <- NA

gta_cellanno <- rbind(gta_cellanno, gta_ref_cellanno)
gta_cellanno <- gta_cellanno[match(colnames(mat_cpkt), rownames(gta_cellanno)), ]
rownames(gta_cellanno) <- colnames(mat_cpkt)

gta_rowrange <- GRanges(cpkt$CNAmat[, 1], IRanges(cpkt$CNAmat[, 2], cpkt$CNAmat[, 2]+1))

gta_mix <- SingleCellExperiment(
  assays=list(logCNA=mat_cpkt, CNA=2^mat_cpkt-1e-3),
  colData=gta_cellanno
)
rowRanges(gta_mix) <- gta_rowrange
gta_mix <- addPCA(gta_mix, npcs=30, slot = 'logCNA') #logCNA for ARTC02, ARTC94
gta_mix <- addUMAP(gta_mix, use_pcs = 1:30) 
gta_mix <- run_phenograph(gta_mix, k=50, use_dims=1:2, reduction='UMAP')
leiden_res <- 0.1
while(TRUE){
  gta_mix <- run_leiden(gta_mix, resolution = leiden_res)
  leiden_str <- sprintf('Cluster_leiden_res_%s', leiden_res)
  leiden_n_cluster <- length(levels(colData(gta_mix)[[leiden_str]]))
  colData(gta_mix)[['icna_clusters']] <- colData(gta_mix)[[leiden_str]]
  if (leiden_n_cluster >= 5) {break()} ## Maybe 1 obs-normal, 1 ref-normal, 2 obs-tumors, 1 obs-tumor/normal
  if (leiden_res >= 2) {break()}  ## Higher resolutin won't help anymore; probably no aneuploidy cells at all.
  colData(gta_mix)[[leiden_str]] <- NULL
  leiden_res <- leiden_res + 0.1
}


p <- drplot(gta_mix, 'UMAP', group.by=leiden_str) +
  scale_color_discrete(
    labels=pretty_table2str(table(colData(gta_mix)[[leiden_str]])))
ggsave2(file.path(dir_copykat, 'copykat_DNA_DR.leiden_cluster'), 
        p, width = 6, height = 5.5)
p <- drplot(gta_mix, 'UMAP', group.by='category', add_label = F)
ggsave2(file.path(dir_copykat, 'copykat_DNA_DR.category'), 
        p, width = 6, height = 5.5)
# p
table(colData(gta_mix)[[leiden_str]], 
      colData(gta_mix)[['copykat_pred_tirosh']])

cli_alert_info('Find clusters without any reference cell; this cluster is probably the cancer cells.')
print(table(colData(gta_mix)[[leiden_str]],
            colData(gta_mix)[['category']]))
tmp <- table(colData(gta_mix)[[leiden_str]],
             colData(gta_mix)[['category']]) %>% as.matrix()
tmp_prop <- table(colData(gta_mix)[[leiden_str]],
                  colData(gta_mix)[['category']]) %>% 
  prop.table(., margin = 1) %>% as.matrix()

print(tmp)
leiden_is_tumor <- tmp[, 'ref'] < 10 & tmp_prop[, 'ref'] <= 1/100
leiden_is_tumor <- ifelse(leiden_is_tumor, 'aneuploid', 'diploid')

copykat_pred_leiden <- colData(gta_mix)[[leiden_str]]
copykat_pred_leiden <- recode_factor(copykat_pred_leiden, !!!leiden_is_tumor)
table(copykat_pred_leiden)
copykat_pred_leiden <- factor(as.character(copykat_pred_leiden), 
                              levels = c('aneuploid', 'diploid'))
colData(gta_mix)[['copykat_pred_leiden']] <- copykat_pred_leiden
idx <- match(Cells(sr3), colnames(gta_mix))
tmp <- as.character(colData(gta_mix)[['copykat_pred_leiden']][idx])
tmp <- replace_na(tmp, 'diploid')
tmp <- factor(as.character(tmp), levels = c('aneuploid', 'diploid'))
sr3[['copykat_pred_leiden']] <- tmp  ## --> to be exported
print(table(sr3[['copykat_pred_leiden']]))

table(colData(gta_mix)[[leiden_str]], 
      colData(gta_mix)[['copykat_pred_leiden']])
p_scatter_dr <- DimPlot(
  sr3, reduction = 'umap', group.by='copykat_pred_leiden', order = 'aneuploid') +
  theme_void() + 
  scale_color_manual(values = color_copykat_ploidy,
                     labels = pretty_table2str(table(sr3$copykat_pred_leiden))) +
  coord_fixed() + labs(title=sample_name) + theme(legend.position = 'bottom')
ggsave(file.path(dir_copykat, 'copykat_pred_leiden-umap.png'), 
       plot = p_scatter_dr, width = 4.2, height = 4)
ggsave(file.path(dir_copykat, 'copykat_pred_leiden-umap.pdf'), 
       plot = p_scatter_dr, width = 4.2, height = 4)  

## Viz icna_leiden clusters
tmp <- colData(gta_mix[, intersect(Cells(sr3), colnames(gta_mix))])[['icna_clusters']]
names(tmp) <- intersect(Cells(sr3), colnames(gta_mix))
tmp <- tmp[match(Cells(sr3), names(tmp))]
sr3$icna_leiden <- tmp
p_scatter_dr <- UMAPPlot(sr3, group.by='icna_leiden') + theme_void(base_size = 20)+
  labs(title=sample_name, color='icna_leiden')+
  coord_equal() + theme(legend.position = 'bottom')
ggsave(file.path(dir_copykat, 'copykat_icna_leiden-umap.png'), 
       plot = p_scatter_dr, width = 4.2, height = 4)
ggsave(file.path(dir_copykat, 'copykat_icna_leiden-umap.pdf'), 
       plot = p_scatter_dr, width = 4.2, height = 4)  
rm(tmp)

#-------------------------- Livnat Prediction --------------------------
cli_rule('Livnat Prediction')
copykat_pred_livnat <- rep('diploid', ncol(sr3))  ## initiate
names(copykat_pred_livnat) <- Cells(sr3)
copykat_pred_livnat <- factor(copykat_pred_livnat, 
                              levels = c('aneuploid', 'diploid'))

copykat_pred_tirosh <- sr3$copykat_pred_tirosh
copykat_pred_leiden <- sr3$copykat_pred_leiden
cells_tumor_cand <- c(as.character(copykat_pred_tirosh) %in% 'aneuploid') | c(as.character(copykat_pred_leiden) %in% 'aneuploid')
cells_tumor_cand[is.na(cells_tumor_cand)] <- FALSE

if (sum(cells_tumor_cand, na.rm = T)==0){
  cat("No tumor cell candidates, so NO need to perform the actual Livnat-prediction.\n")
} else {
  ## Proceed the Livnat-prediction
  cells_tumor_cand <- Cells(sr3)[cells_tumor_cand]
  
  cat(length(cells_tumor_cand), '/', ncol(sr3), 'cells are malignant candidates.\n')

  library(SummarizedExperiment)
  library(future.apply)
  lib_tnbc_tumor <- readRDS("data/lib/tcga/brca/tnbc_primary_tumor.se.rds")
  lib_tnbc_normal <- readRDS("data/lib/tcga/brca/tnbc_solid_normal.se.rds")
  
  mat_sc <- GetAssayData(sr3[, cells_tumor_cand], assay='RNA', slot='data')
  features_shared <- VariableFeatures(sr3)
  
  features_shared <- intersect(
    features_shared, 
    intersect(rowData(lib_tnbc_tumor)[['external_gene_name']],
              rowData(lib_tnbc_normal)[['external_gene_name']]))
  
  str(features_shared)
  mat_sc <- mat_sc[features_shared, ]

  adhoc_prepare_lib_mat <- function(se, features, feature_type='Symbol', 
                                    assay='logFPKM'){
    data <- assays(se)[[assay]]
    features_se <- rowData(se)[['external_gene_name']] ## symbol
    features_se <- make.unique(features_se)
    features <- intersect(features, features_se)
    idx <- match(features, features_se)
    data <- data[idx, ]
    rownames(data) <- features
    return(data)
  }
  mat_tnbc_tumor <- adhoc_prepare_lib_mat(lib_tnbc_tumor, features_shared)
  mat_tnbc_normal <- adhoc_prepare_lib_mat(lib_tnbc_normal, features_shared)
  ## Livnat correlation-based prediction
  adhoc_wilcox_a_vs_b_is_bigger <- function(a, b){
    obj <- wilcox.test(a, b, alternative = 'greater')
    if (obj$p.value < 0.1) {return(TRUE)}
    return(FALSE)
  }
  df_livnat <- future_apply(mat_sc, 2, function(v) {
    cor_over_tumor_i <- as.numeric(cor(mat_tnbc_tumor, as.numeric(v)))
    cor_over_normal_i <- as.numeric(cor(mat_tnbc_normal, as.numeric(v)))
    pval_i <- wilcox.test(cor_over_normal_i, cor_over_tumor_i, 'g')
    pval_i <- pval_i$p.value
    # boxplot(list('norm'=cor_over_normal_i, 'tumor'=cor_over_tumor_i))
    return(c(`avgcor_tumor_tissue`=mean(cor_over_tumor_i, rm.na=T), 
             `avgcor_normal_tissue`=mean(cor_over_normal_i, rm.na=T),
             `pval` = pval_i))
  })
  is_cand_normal <- df_livnat['pval', ] < 0.05
  print(table(is_cand_normal))
  ## todo: two hists to see the corr
  is_cand_normal[is.na(is_cand_normal)] <- TRUE
  cells_tumor_confirmed <- cells_tumor_cand[!is_cand_normal]
  cat(sum(is_cand_normal), 'tumor candidateds are trimmed due to their similarity to bulk normal tissue.\n')
  cat(length(cells_tumor_confirmed), 'tumor candidated are further confirmed.\n')
  
  if (sum(!is_cand_normal) > 0) {
    copykat_pred_livnat[cells_tumor_confirmed] <- 'aneuploid'
  }
  print(table(copykat_pred_livnat))
}


sr3[['copykat_pred_livnat']] <- copykat_pred_livnat  ## to-be exported
print(table(sr3[['copykat_pred_livnat']]))


p_scatter_dr <- DimPlot(
  sr3, reduction = 'umap', group.by='copykat_pred_livnat', order = 'aneuploid') +
  theme_void() + 
  scale_color_manual(values = color_copykat_ploidy,
                     labels = pretty_table2str(table(sr3$copykat_pred_livnat))) +
  coord_fixed() + labs(title=sample_name) + theme(legend.position = 'bottom')

ggsave(file.path(dir_copykat, 'copykat_pred_livnat-umap.pdf'), 
       plot = p_scatter_dr + theme_pubr(), width = 4.2, height = 4)  

idx <- match(colnames(gta_mix), Cells(sr3))
tmp <- sr3$copykat_pred_livnat[idx]
tmp[is.na(tmp)] <- 'diploid'
colData(gta_mix)[['copykat_pred_livnat']] <- tmp
rm(idx); rm(tmp)



#-------------------------- Export Results --------------------------
# copykat_pred_default ## Ruli's default hclust=2
# copykat_pred_tirosh  ## CNV score v.s. CNV corr
# copykat_pred_leiden  ## clustering on iCNAs
# copykat_pred_livnat  ## corr with tumor/normal bulk RNA samples

final_report <- FetchData(
  sr3, c('copykat_pred_default', 'copykat_pred_tirosh', 
         'copykat_pred_leiden', 'copykat_pred_livnat'))
write_tsv(rownames_to_column(final_report, 'cellnames'), 
          file.path(dir_copykat, 'copykat_pred_report.tsv'))
write_rds(final_report, file.path(dir_copykat, 'copykat_pred_report.seurat3_meta.rds'))
write_rds(gta_mix, file.path(dir_copykat, 'copykat_mix.sce.rds'))
cli_rule('[DONE]')
timestamp()


#-------------------------- Final Visualization --------------------------  
cli_rule('Final Heatmap Viz')
set.seed(22)
hm_obj_sc <- plot_heatmap_sc_manual(
  gta_mix, clip=c(-1, 1),
  cell_group_by = leiden_str,
  anno_rows_category = c('category',
                         'copykat_pred_default', 
                         'copykat_pred_tirosh',
                         'copykat_pred_leiden',
                         'copykat_pred_livnat'))
pdf(file.path(dir_copykat, 'copykat_heatmap3.pdf'), width = 13, height = 10)
draw(hm_obj_sc,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

