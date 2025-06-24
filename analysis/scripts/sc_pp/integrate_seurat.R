#---------------------------
# Perform the default seurat intergration              ----    
#---------------------------
library(Seurat)
library(patchwork)
library(tidyverse); library(readr)

#------ Inputs ------

cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
  fpath_sr3 <- cmdargs[1]
  n_pc_nn <- cmdargs[2]
} else {
  fpath_sr3 <- ''
  n_pc_nn <- Inf
}


#------------------ ~~~ Begin ~~~ --------------------

dir_res <- file.path(dirname(fpath_sr3), sprintf('integrate_seurat_%spc', n_pc_nn))
fs::dir_create(dir_res)
topath_sr3 <- file.path(dir_res, 'ready.sr3.rds')
topath_metadata <- ''

message('------ Read object')
sr3 <- read_rds(fpath_sr3)

message('------ Preprocessing')
DefaultAssay(sr3) <- 'RNA'

message('------ Prepare list of objects')
if (F) {
  # correct by patient
  sr3_list <- SplitObject(sr3, split.by = "patient")
  k.filter <- 50
  k.weight <- 30
}
if (F) {
  # by random batch
  n_random_batch <- 10
  sr3$random_batch <- factor(sample(1:n_random_batch, size=ncol(sr3), replace = T))
  sr3_list <- SplitObject(sr3, split.by = "random_batch")
  k.filter <- 200
  k.weight <- 100
}
if (T) {
  # correct patient group
  n_random_batch <- 5
  sr3$random_batch <- factor(
    c(as.numeric(as.factor(sr3$patient)) %% n_random_batch) + 1
  ) 
  sr3_list <- SplitObject(sr3, split.by = "random_batch")
  k.filter <- 100
  k.weight <- 100
}


sr3_list <- lapply(sr3_list, function(x){ 
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
})

message('------ SelectIntegrationFeatures')
features <- SelectIntegrationFeatures(object.list = sr3_list)

message('------ SelectIntegrationFeatures')
sr3_anchors <- FindIntegrationAnchors(
  object.list = sr3_list, anchor.features = features, k.filter = k.filter)
write_rds(sr3_anchors, paste0(topath_sr3, '.anchor'))

sr3_int <- IntegrateData(anchorset = sr3_anchors, k.weight = k.weight)

DefaultAssay(sr3_int) <- "integrated"
write_rds(sr3_int, topath_sr3)

message('------ Perform an integrated analysis')
sr3_int <- ScaleData(sr3_int, verbose = FALSE)
sr3_int <- RunPCA(sr3_int, npcs = 300, verbose = FALSE)
sr3_int <- RunUMAP(sr3_int, reduction = "pca", dims = 1:n_pc_nn)
sr3_int <- FindNeighbors(sr3_int, reduction = "pca", dims = 1:n_pc_nn)
sr3_int <- FindClusters(sr3_int, resolution = 0.8)
write_rds(sr3_int, topath_sr3)
write_rds(x = sr3_int@meta.data, path = file.path(dir_res, 'sr3_metadata.df.rds'))

library(ggpubr)
# sr3_int <- read_rds(topath_sr3)
p <- UMAPPlot(sr3_int, group.by ='patient') + rremove('legend') + theme(aspect.ratio = 1)
ggsave(file.path(dir_res, 'dr.patient.pdf'), p, width = 7, height = 7)


cat('\n==DONE==\n')