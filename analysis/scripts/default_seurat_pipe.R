#--------------------------
# Prepare analysis-ready seurat object (mostly use RNA assya)
#--------------------------

suppressPackageStartupMessages({
  library(Seurat); library(S4Vectors)
  library(tidyverse); library(ggplot2); library(ggpubr); library(patchwork)
  theme_set(theme_pubr(base_size = 11, legend='right')); 
  library(cli); library(tictoc); library(glue); library(scales); library(tools)
  library(purrr); library(stringr)
  library(clustree)
  library(future)
  options(future.globals.maxSize = 10 *1024^3)
  plan("multisession", workers = 10)
  plan("sequential")
  library(RhpcBLASctl)
  RhpcBLASctl::blas_set_num_threads(10)
})

args <- commandArgs(trailingOnly = TRUE)
do_scale_nUMI <- NULL
CCC_plan <- NULL
n_pc <- NULL
n_pc_nn <- NULL
assay_use <- NULL
n_var_genes <- NULL


#-------------------------- START --------------------------    
if (!is_empty(args)) {
  f_sr3 <- args[[1]]
  dir_res <- args[[2]]
  CCC_plan <- try(args[[3]]) # No=0,  std=1,  alt=2
  n_pc <- try(as.numeric( args[[4]] ))
  n_pc_nn <- try(as.numeric( args[[5]] ))
  n_var_genes <- try(as.numeric(args[[6]]))
  assay_use <- try(args[[7]])
}

if ('try-error' %in% class(do_scale_nUMI) | is.null(do_scale_nUMI)) {do_scale_nUMI <- TRUE}
if ('try-error' %in% class(CCC_plan) | is.null(CCC_plan)) {CCC_plan <- 1} else {CCC_plan <- as.numeric(CCC_plan)}
if ('try-error' %in% class(n_pc) | is.null(n_pc)) { n_pc <- 300}
if ('try-error' %in% class(n_pc_nn) | is.null(n_pc_nn)) { n_pc_nn <- 200}
if ('try-error' %in% class(assay_use) | is.null(assay_use)) { assay_use <- 'RNA'}
if ('try-error' %in% class(n_var_genes) | is.null(n_var_genes)) { n_var_genes <- 2000 }


stopifnot(file.exists(f_sr3))

f_sr_out <- file.path(dir_res, 'ready.sr3.rds')
f_cellmeta <- file.path(dir_res, 'sr3_metadata.df.rds')

fs::dir_create(dir_res)
#--------------------------
# Make analysis-ready
#--------------------------
do_CCC_alt <- CCC_plan==2
do_CCC <- CCC_plan==1


stopifnot(!all(c(do_CCC_alt, do_CCC)))


if (FALSE){
  hvg_blacklist <- read_rds(
    file.path(
      '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat40',
      'lv01.aneuploidy_tri_type.aneuploid.10xV3.RNA.DeDb', 'non_patient_specific_HVGs',
      'patient_private_genes.vector.rds'))
  hvg_blacklist <- read_rds(
    file.path(
      '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat40',
      'lv01.aneuploidy_tri_type.aneuploid.10xV3.RNA.DeDb', 'non_patient_specific_HVGs',
      'patient_specific_genes.vector.rds'))
}
hvg_use <- NULL
hvg_blacklist <- NULL

#------ Cell types genes ------
if (FALSE){
  lib_celltype_rna_markers <- "/volumes/USR1/yyan/project/tumor_plasticity/lib/cellType_markerGenes-all_colored-NoDup-v6.csv"
  lib_celltype_rna_markers <- suppressMessages(read_csv(
    lib_celltype_rna_markers, col_names = T))
  manual_celltypes_color <- deframe(lib_celltype_rna_markers[, c('Type', 'Color')])
  level_use_min = 10; level_use_max = 30
  lib_celltypes_use <- lib_celltype_rna_markers %>%
    dplyr::filter(Level>=level_use_min, Level<=level_use_max, Type != 'Unknown', Type != 'Dead')
  adhoc_genes_str_to_vector <- function(s){
    s <- stringr::str_remove_all(s, ' ')
    sapply(s, function(ss) unique(unlist(strsplit(ss, split=','))), simplify = T)
  }
  markers_use <- adhoc_genes_str_to_vector(lib_celltypes_use$Genes)
  names(markers_use) <- lib_celltypes_use$Type
  hvg_use <- unname(unlist(markers_use))
  str(hvg_use)
}

    
#------ Run Analysis ------
  
if (! file.exists(f_sr_out)){
  sr3 <- read_rds(f_sr3)
  DefaultAssay(sr3) <- assay_use
  sr3 <- DietSeurat(object = sr3, assays = assay_use)
  print(sr3)
  #------ Find HVG ------
  if (is_empty(hvg_use)){
    sr3 <- FindVariableFeatures(sr3, nfeatures = n_var_genes)
  } else {
    cat('Manual defined hvg.')
    VariableFeatures(sr3) <- hvg_use
  }
  
  if (!is_empty(hvg_blacklist)){
    tmp <- VariableFeatures(sr3)
    cat(sum(tmp %in% hvg_blacklist), 'genes are patient private.\n')
    tmp <- setdiff(tmp, hvg_blacklist)
    VariableFeatures(sr3) <- tmp
  }
  print(sr3)
  #------ Cell Scoring ------
  if (!'Phase' %in% colnames(sr3[[]])){
    sr3 <- CellCycleScoring(
      sr3, s.features = Seurat::cc.genes.updated.2019$s.genes,
      g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
      assay='RNA',
      set.ident = FALSE)
  }
  #------ Scale - cellcycle or not ------
  vars.to.regress.str <- NULL
  if (do_CCC_alt) {
    cat('Regress out CC.Difference')
    #------ CCC_alt Keep cycling v.s. non-cycling ------
    sr3$CC.Difference <- sr3$S.Score - sr3$G2M.Score
    vars.to.regress.str <- "CC.Difference"
  } 
  if (do_CCC){
    cat('Standard Cell Cycle Correction')
    #------ CCC ------
    vars.to.regress.str <- c("S.Score", "G2M.Score")
  }; cat('\n')
  if (!do_CCC & !do_CCC_alt){
    #------ Simple Re-Scale ------
    vars.to.regress.str <- NULL
  }
  if (do_scale_nUMI) {
    vars.to.regress.str <- c(vars.to.regress.str, 'nCount_RNA')  
  }
  str(VariableFeatures(sr3))
  sr3 <- ScaleData(sr3, assay = assay_use, 
                  features = VariableFeatures(sr3),
                  vars.to.regress = vars.to.regress.str)
  #------ Dimention Reduction ------
  n_pc <- min(ncol(sr3)-1, length(VariableFeatures(sr3))-1, n_pc)
  n_pc_nn <- min(ncol(sr3)-1, length(VariableFeatures(sr3))-1, n_pc_nn)
  
  sr3 <- RunPCA(sr3, npcs = n_pc, assay = assay_use, features = VariableFeatures(sr3))
  p_elbow <- ElbowPlot(sr3, ndims=n_pc) + geom_vline(xintercept = n_pc_nn)
  ggsave2(file.path(dir_res, 'pca_elbow'), p_elbow, width = 6, height = 6)
  if (ncol(sr3) > 50) {
    sr3 <- RunUMAP(sr3, reduction = 'pca', assay = assay_use,
                  dims = 1:n_pc_nn, spread = 1.2)
  }

  if (F){
    sr3 <- RunICA(object = sr3, nics=n_pc, assay = assay_use)
  }  

  #------ Export sr3 ------
  write_rds(sr3, f_sr_out)
  df_cellmeta <- sr3[[]]
  write_rds(df_cellmeta, f_cellmeta)
} else {
  message('Just read the existing object. Skipped all the pre-processing. ')
  sr3 <- read_rds(f_sr_out)
  df_cellmeta <- read_rds(f_cellmeta)
  
}

#------ Build NN Graph ------
sr3 <- FindNeighbors(
  sr3, dims = 1:n_pc_nn, assay = assay_use,
  reduction = 'pca')

#------ Clustering ------
## clear previous snn result
idx <- grepl(pattern = '_snn_res', x=colnames(sr3[[]]))
for (i in colnames(sr3[[]])[idx]) {
  sr3[[i]] <- NULL
}
snn_res_max <- 1
if (ncol(sr3) < 100) {snn_res_max <- 1}
for (snn_res_i in seq(from=0.8, to=snn_res_max, by = .2)) {
  cat(snn_res_i, '... ')
  try( sr3 <- FindClusters(sr3, resolution = snn_res_i, verbose = F) )
  snn_str_i <- sprintf('%s_snn_res.%s', assay_use, snn_res_i)
  if (snn_str_i %in% colnames(sr3[[]])) {
    sr3[[snn_str_i]] <- Idents(sr3) ## The cluster order is human-readable
  }
}; cat('\n')

#------ Export sr3 ------
write_rds(sr3, f_sr_out)
df_cellmeta <- cbind(
  FetchData(sr3, c('UMAP_1', 'UMAP_2', 'tSNE_1', 'tSNE_2', 'PC_1', 'PC_2')),
  sr3[[]])
write_rds(df_cellmeta, f_cellmeta)

#------ Evaluate clustering resolution ------
library(clustree)
p <- clustree(sr3, prefix = sprintf("%s_snn_res.", assay_use))
ggsave2(file.path(dir_res, 'clustree'), p, width = 10, height = 15)
} else {
  sr3 <- read_rds(f_sr_out)
  df_cellmeta <- read_rds(f_cellmeta)
  
}

cat('DONE!\n')
timestamp()
