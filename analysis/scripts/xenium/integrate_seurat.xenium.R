#---------------------------
# Perform the default seurat intergration              ----    
#---------------------------
library(Seurat)
library(patchwork)
library(tidyverse); library(readr)
library(ggpubr)
source('~/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R')
source('~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R')
my_scatter_themevoid <- theme_pubr(base_size = 6, legend = 'right') %+replace% theme(
  aspect.ratio = 1, 
  axis.text=element_blank(), 
  axis.title=element_blank(), 
  axis.ticks=element_blank(), 
  panel.border = element_rect(fill = NA, linewidth=rel(.5)), 
  axis.line = element_blank())
#------ fibroblast ------
# fpath_sr3 <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Fibro/ready.sr3.rds'
# n_pc_nn <- 50
# #------ stromal compartment ------
# fpath_sr3 <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/TNBC_stromal_cells/ready.sr3.rds'
# n_pc_nn <- 300
# 
# #------ atlas cells (1e3 cells per state) ------
# fpath_sr3 <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/downto_1000_by_cell_state_paper/ready.sr3.rds'
# n_pc_nn <- 100

fpath_sr3 <- '/volumes/USR1/yyan/project/tnbc_xenium/data_merged_xenium5kPlus/celltype_T/ready.seurat.rds'
n_pc_nn <- 50

cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {
  fpath_sr3 <- cmdargs[1]
  n_pc_nn <- as.numeric(cmdargs[2])
}
#------------------ ~~~ Begin ~~~ --------------------

dir_res <- file.path(dirname(fpath_sr3), sprintf('integrate_seurat_%spc', n_pc_nn))
fs::dir_create(dir_res)
topath_sr3 <- file.path(dir_res, 'ready.sr3.rds')
topath_metadata <- file.path(dir_res, 'sr3_metadata.df.rds')

# if (file.exists(topath_sr3)) {
#   cat('No need to run')
#   q(save = 'no')
# }

sr3 <- read_rds(fpath_sr3)
if ('sample' %in% colnames(sr3@meta.data) ) {
  pal_sample <- init_pal_d(as.character(unique(sr3$sample)))
  if ('umap' %in% names(sr3@reductions)) {
    p <- DimPlot(sr3, group.by ='sample', reduction = 'umap', raster = T,
                 label = TRUE, cols = pal_sample, 
                 raster.dpi = c(1024, 1024), pt.size = 3, shuffle = T) +
      my_scatter_themevoid
    ggsave(file.path(dir_res, 'dr.before_integration.sample.pdf'),
           p + rremove('legend') , width = 7, height = 7)
    ggsave(file.path(dir_res, 'legend.sample.pdf'),
           as_ggplot(get_legend(p)), width = 3, height = 7)
  }
}

if (! file.exists(topath_sr3) ) {
# if (T) {
  message('------ Read object')
  
  message('------ Preprocessing')
  DefaultAssay(sr3) <- 'Xenium'
  
  message('------ Prepare list of objects')
  if (F) {
    # correct sample
    sr3_list <- SplitObject(sr3, split.by = "sample")
    k.filter <- NA
    k.weight <- 100
  }
  if (F) {
    # simply random batch -- not recommend
    n_random_batch <- 10
    sr3$random_batch <- factor(sample(1:n_random_batch, size=ncol(sr3), replace = T))
    sr3_list <- SplitObject(sr3, split.by = "random_batch")
    k.filter <- NA
    k.weight <- 100
  }
  if (T) {
    # correct sample group -- each sample is well distributed to all batches
    n_random_batch <- 5
    cell2cat <- structure(as.character(sr3$sample), 
                          names=Cells(sr3))
    set.seed(1026); cell2cat <- cell2cat[sample(seq_len(length(cell2cat)), 
                                                size=length(cell2cat))]
    random_batch <- c(as.numeric(as.factor(cell2cat)) %% n_random_batch) + 1
    names(random_batch) <- names(cell2cat)
    
    random_batch <- random_batch[Cells(sr3)]
    sr3$random_batch <- as.factor(random_batch)
    
    sr3_list <- SplitObject(sr3, split.by = "random_batch")
    k.filter <- NA
    k.weight <- 100
  }
  
  
  sr3_list <- lapply(sr3_list, function(x){ 
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nrow(x))
    # VariableFeatures(x) <- rownames(x)
    x
  })
  num_hvg_list <- sapply(sr3_list, function(x) length(VariableFeatures(x)))
  # n_integration_features <- min(c(5000, num_hvg_list))
  n_integration_features <- min(num_hvg_list)
  cat('use ', n_integration_features, ' featuers to integration.\n')
  message('------ SelectIntegrationFeatures')
  features <- SelectIntegrationFeatures(
    object.list = sr3_list, 
    nfeatures = n_integration_features) ## 5k gene panel
  
  # fvf.nfeatures
  
  message('------ FindIntegrationAnchors')
  if (! file.exists(paste0(topath_sr3, '.anchor'))) {
    sr3_anchors <- FindIntegrationAnchors(
      object.list = sr3_list, 
      anchor.features = features, 
      k.filter = k.filter, 
      reduction  = "cca")
    write_rds(sr3_anchors, paste0(topath_sr3, '.anchor'))
  } else {
    sr3_anchors <- read_rds(paste0(topath_sr3, '.anchor'))
  }
  
  message('------ IntegrateData')
  sr3_int <- IntegrateData(anchorset = sr3_anchors, k.weight = k.weight)
  
  DefaultAssay(sr3_int) <- "integrated"
  write_rds(sr3_int, topath_sr3)
} else {
  sr3_int <- readRDS(topath_sr3)
}  
message('------ Perform an integrated analysis')
print(dim(sr3_int))
if (!file.exists(file.path(dir_res, 'feature_names_integration.rds'))) {
  sr3_int <- ScaleData(sr3_int, verbose = FALSE)
  sr3_int <- RunPCA(sr3_int, npcs = n_pc_nn, verbose = FALSE)
  sr3_int <- RunUMAP(sr3_int, reduction = "pca", dims = 1:n_pc_nn)
  sr3_int <- FindNeighbors(sr3_int, reduction = "pca", dims = 1:n_pc_nn)
  sr3_int <- FindClusters(sr3_int, resolution = 0.8)
  write_rds(sr3_int, topath_sr3)
  write_rds(x = sr3_int@meta.data, path = topath_metadata)
  write_rds(rownames(sr3_int), file.path(dir_res, 'feature_names_integration.rds'))
  write_lines(rownames(sr3_int), file.path(dir_res, 'feature_names_integration.txt'))
} else {
  sr3_int <- readRDS(topath_sr3)
}

#------ clustering ------
idx <- grepl(pattern = '_snn_res', x=colnames(sr3_int[[]]))
for (i in colnames(sr3_int[[]])[idx]) {sr3_int[[i]] <- NULL}
snn_res_max <- 1
if (ncol(sr3_int) < 100) {snn_res_max <- .2}

for (snn_res_i in seq(from=0.2, to=snn_res_max, by = .2)) {
  cat(snn_res_i, '... ')
  try( sr3_int <- FindClusters(
    sr3_int, algorithm = 3, 
    resolution = snn_res_i, verbose = F) )
  snn_str_i <- sprintf('%s_snn_res.%s', DefaultAssay(sr3_int), snn_res_i)
  if (snn_str_i %in% colnames(sr3_int[[]])) {
    sr3_int[[snn_str_i]] <- Idents(sr3_int) ## The cluster order is human-readable
  }
}; cat('\n')
write_rds(sr3_int, topath_sr3)
write_rds(x = sr3_int@meta.data, path = topath_metadata)

#------------------ ~~~ Viz ~~~ --------------------

#------ viz technical features ------

for (feature_i in c('nCount_Xenium', 'nFeature_Xenium',
                    'nCount_RNA', 'nFeature_RNA',
                    'nCount_SCT', 'nFeature_SCT', 
                    'EPCAM')) {
  if (! (feature_i %in% colnames(sr3_int@meta.data) | feature_i %in% rownames(sr3_int))) {
    cat('requested feature ', feature_i, 'is not available...\n')
    next()
  }
  message(feature_i)
  
  dir_snippet_viz <- dir_res
  pal_v <- c('lightgrey', 'blue')## either c('low_color', 'high_color') or just NULL (to be change outside)
  viz_which_val <- feature_i
  
  pa <- FeaturePlot(
    sr3_int, features = viz_which_val, 
    reduction = 'umap',
    max.cutoff = 'q95', min.cutoff = 'q5', 
    raster = T, pt.size = 3, raster.dpi = c(1024, 1024),
    cols = pal_v) + 
    labs(caption = sprintf('%s cells', ncol(sr3_int))) + 
    my_scatter_themevoid 
  ruok::ggsave2(file.path(dir_snippet_viz, 
                          sprintf('featureplot.%s', viz_which_val)), 
                pa + rremove('legend'), width=7, height = 7)
  ggsave(file.path(dir_snippet_viz, 
                   sprintf('legend.featureplot.%s.pdf', viz_which_val)), 
         as_ggplot(get_legend(pa)), width=1.5, height = 3, useDingbats = F)
  
}



#------ viz ident ------

snn_res_i <- 0.2
snn_str_i <- sprintf('%s_snn_res.%s', DefaultAssay(sr3_int), snn_res_i)
snn_str_opts <- colnames(sr3_int@meta.data)[grepl(sprintf('%s_snn_res', DefaultAssay(sr3_int)), colnames(sr3_int@meta.data))]
viz_z_opts <- c(snn_str_opts, c('sample', 'patient'))
viz_z_opts <- intersect(viz_z_opts, colnames(sr3_int@meta.data))
viz_z_opts
for (snn_str_i in viz_z_opts) {
  cat(snn_str_i, '... ')
  sr3_int$seurat_clusters <- sr3_int@meta.data[, snn_str_i]
  
  if (!'factor' %in% class(sr3_int$seurat_clusters)) {sr3_int$seurat_clusters <- as.factor(sr3_int$seurat_clusters)}
  
  pal_snn <- structure(
    Seurat::DiscretePalette(n=length(levels(sr3_int$seurat_clusters)), palette = 'parade'),
    names=levels(sr3_int$seurat_clusters))
  
  dir_snippet_viz <- dir_res
  viz_what <- snn_str_i
  pal_z <- pal_snn
  pa <- DimPlot(
    sr3_int, reduction = 'umap', label = T, 
    raster = T, pt.size = 3, raster.dpi = c(1024, 1024),
    group.by = viz_what, shuffle = T, cols = pal_z) + 
    labs(caption = sprintf('%s cells', ncol(sr3_int))) + 
    my_scatter_themevoid
  ruok::ggsave2(file.path(dir_res, 
                          sprintf('dimplot.%s', viz_what)), 
                pa + rremove('legend'), width=7, height = 7)
  ggsave(file.path(dir_snippet_viz, 
                   sprintf('legend.dimplot.%s.pdf', viz_what)), 
         as_ggplot(get_legend(pa)), width=1.5, height = 3, useDingbats = F)

}; cat('[done]\n')


# library(ggpubr)
# sr3_int <- read_rds(topath_sr3)
# p <- DimPlot(sr3_int, group.by ='sample', reduction = 'umap', raster = T, 
#              raster.dpi = c(1024, 1024), pt.size = 3, shuffle = T) + 
#   my_scatter_themevoid
# ggsave(file.path(dir_res, 'dr.sample.pdf'), 
#        p + rremove('legend') , width = 7, height = 7)
# ggsave(file.path(dir_res, 'legend.sample.pdf'), 
#        as_ggplot(get_legend(p)), width = 3, height = 7)
# p <- UMAPPlot(sr3_int, group.by ='cell_state_paper', label=T) + theme(aspect.ratio = 1)
# ggsave(file.path(dir_res, 'dr.cell_state_paper.pdf'), p, width = 7, height = 7)
# p <- UMAPPlot(sr3_int, group.by ='ecotrait_feature', label=T) + theme(aspect.ratio = 1)
# ggsave(file.path(dir_res, 'dr.ecotrait_feature.pdf'), p, width = 7, height = 7)


cat('\n==DONE==\n')
