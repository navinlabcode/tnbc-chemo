#---------------------------
# Visualize features (genes / signatures) across samples              ----    
# Yun Yan
#---------------------------
library(Seurat)
library(colorspace); library(ggplot2); library(scales); library(ggpubr)
library(readr)
library(ggbeeswarm)

theme_set(theme_pubr())
setwd('/volumes/USR1/yyan/project/tnbc_spatial')
source('util.visium.R')
source('util.ruok.R')

dir_res <- './result/viz_features_across_samples'
fs::dir_create(dir_res)

run_st_list <- c('ART23', 'ART10')#, 'ART43') [ARTC10=ARC2, ARTC23=ARC3]
# run_st_list <- c('ART23', 'ART65')#, 'ART43') [ARTC10=ARC2, ARTC23=ARC3]
# run_st_list <- c('ART23', 'ART43')#, 'ART43') [ARTC10=ARC2, ARTC23=ARC3]

features_list <- c('MTDH', 'CDKN2A') #[MTDH=ARC2, CDKN2A=ARC3]
features_list <- c('MTDH', 'CDKN2A', 'CKS1B', 'MYC', 'GADD45A', 'CCND1', 'GATA3')
features_list <- c('HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DMA',
                   'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 
                   'TOP2A', 'MKI67', 'CENPF', 
                   'PCNA', 'TYMS', 'CLSPN', 
                   'ISG15', 'IFI44', 'IFIT3', 
                   'HSPA5', 'DNAJB11', 'SERP1', 
                   'TAGLN', 'ACTA2', 'MYLK', 
                   'S100A6', 'S100A10', 'LGALS3', 
                   'TMSB4X', 'ANXA1', 'KRT19', 'KRT81', 
                   'HIST1H4C', 'PCLAF', 'ATAD2', 'DUT', 'TK1', 'H2AFZ', 
                   'UBE2T', 'TUBB', 'RRM2', 'RAD51AP1', 'DEK'
                   )
# features_list <- c('HSPA5', 'DNAJB11', 'SERP1', 
#                    'TAGLN', 'ACTA2', 'MYLK', 
#                    'S100A6', 'S100A10', 'LGALS3')

f_vsm_list <- file.path('./rds', run_st_list, 'vsm.rds')
names(f_vsm_list) <- run_st_list; all(file.exists(f_vsm_list))

#------ Read ------
## all spots
vsm_list <- lapply(f_vsm_list, read_rds)
## only selected spots
z <- 'confidence'; z_use <- 'High_purity'
lapply(vsm_list, function(o) {
  table(o[[z]])
})
vsm_s_list <- lapply(vsm_list, function(o) {
  Idents(o) <- z
  if (sum(Idents(o) == z_use)==0){
    return(subset(o, idents = 'Medium_purity'))
  }
  return(subset(o, idents = z_use))
})

#------ Parameters ------
  
ASSAY_USE <- 'Spatial'
# ASSAY_USE <- 'SCT'
n_run <- length(vsm_s_list)
covv <-  c('nCount_Spatial', 'nFeature_Spatial', 'nCount_SCT', 'nFeature_SCT')

#------ Output folders ------

for (j in seq_along(features_list)) {
  fs::dir_create(file.path(dir_res, features_list[j]))
}
for (x in covv) {
  fs::dir_create(file.path(dir_res, x))
}

#------ Possible confounding fators such as nUMIs ------

for (covv_use in covv) {
  fs::dir_create(dir_res, covv_use)
  covv_value <- sapply(vsm_s_list, function(o) {
    o@meta.data[, covv_use]
  })
  covv_value_df <- enframe_list(covv_value, name = 'run', value=covv_use)
  p <- ggplot(covv_value_df, aes_string(x='run', y=covv_use)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=2, shape=1) + 
    stat_compare_means()
  ggsave(file.path(dir_res, covv_use, sprintf('vlnplot.%s.%s.%s.pdf', covv_use, z, paste(z_use, collapse = '_'))), 
         p, height = 4, width = 2*n_run, useDingbats = F)  
}
features_list[!features_list %in% rownames(vsm_list[[1]])]
#------ See if gene expression is different across samples ------
for (j in seq_along(features_list)) {
  g_use <- features_list[j]
  message(g_use)
  g_value <- sapply(vsm_s_list, function(o) {
    o@assays[[ASSAY_USE]]@data[g_use, ]
  })
  g_val_df <- enframe_list(g_value, name = 'run', value=g_use)
  pretty_table2str(table(g_val_df$run))
  head(g_val_df)
  colnames(g_val_df) <- make.names(colnames(g_val_df))
  p <- ggplot(g_val_df, aes_string(x='run', y=make.names(g_use))) + 
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=2, shape=16, color='peru', alpha=0.5) + 
    stat_compare_means() + 
    scale_x_discrete(labels=pretty_table2str(table(g_val_df$run))) +
    theme_pubr(base_size = 6) +
    theme(axis.ticks.length=unit(0.06, "inch"))
  print(p)
  Sys.sleep(1)
  ggsave(file.path(dir_res, g_use, sprintf('vlnplot.%s.across_runs.pdf', ASSAY_USE)), 
         p, height = 3, width = 1 * n_run, useDingbats = F)
}
#------ Get gene expression value range ------
global_feature_value_min <- rep(-Inf, length(features_list))
global_feature_value_max <- rep(Inf, length(features_list))
names(global_feature_value_max) <- names(global_feature_value_min) <- features_list
for (j in seq_along(features_list)) {
  g_use <- features_list[j]
  message(g_use)
  g_value <- sapply(vsm_s_list, function(o) {
    o@assays[[ASSAY_USE]]@data[g_use, ]
  })
  g_val_df <- enframe_list(g_value, name = 'run', value=g_use)
  global_feature_value_min[g_use] <- min(g_val_df[[g_use]])
  global_feature_value_max[g_use] <- quantile(g_val_df[[g_use]], 0.5)
}
print(global_feature_value_min)
print(global_feature_value_max)
#------ Viz high/med/low tumor purity spots ------
  
z <- 'confidence'
fs::dir_create(file.path(dir_res, z))
for (i in seq_along(f_vsm_list)) {
  # vsm <- read_rds(f_vsm_list[i]); print(vsm)
  vsm <- vsm_list[[i]]; print(vsm)
  p1 <- SpatialDimPlot(vsm, group.by=z, pt.size.factor=3) +
    scale_fill_viridis_d(direction = -1) +
    theme(aspect.ratio = get_visium_asp_ratio(vsm), legend.position = 'right')
  ggsave2(
    file.path(dir_res, z, sprintf('SpatialDimPlot.%s.%s', z, names(f_vsm_list)[i])),
    p1, width = 5+1, height = 5 * get_visium_asp_ratio(vsm))
  
  vsm_s <- vsm_s_list[[i]]; print(vsm_s)
  p1 <- SpatialDimPlot(vsm_s, group.by=z, pt.size.factor=3, alpha=0.5) + 
    scale_fill_viridis_d(direction = -1) + 
    theme(aspect.ratio = get_visium_asp_ratio(vsm_s), legend.position = 'right')
  ggsave2(
    file.path(dir_res, z, sprintf('SpatialDimPlot.%s_selected.%s', z, names(f_vsm_list)[i])),
    p1, width = 5+1, height = 5 * get_visium_asp_ratio(vsm_s))
}

#------ Viz genes ------

## all spots
if (F) { ## deprecated -- no big use
  for (i in seq_along(f_vsm_list)) {
    vsm <- vsm_list[[i]]; message(names(f_vsm_list)[i])
    for (j in seq_along(features_list)) { ## this is to reduce num of for loops
      message(features_list[j])
      p2 <- SpatialFeaturePlot(vsm, features = features_list[j], slot='data') +
        theme(aspect.ratio = get_visium_asp_ratio(vsm), legend.position = 'right') +
        scale_fill_distiller(palette='Spectral')
      ggsave2(
        file.path(dir_res, features_list[j], 
                  sprintf('SpatialFeaturePlot.%s.%s.all_spots', names(f_vsm_list)[i], features_list[j])),
        p2, width = 5+1, height = 5 * get_visium_asp_ratio(vsm))
    }
  }
}
# Idents(vsm) <- z; table(Idents(vsm))
# z_use <- 'High_purity'
# vsm_s <- subset(vsm, ident=z_use, features = features_list)

## selected spots

for (i in seq_along(f_vsm_list)) {
  vsm_s <- vsm_s_list[[i]]
  message(names(f_vsm_list)[i])
  DefaultAssay(vsm_s) <- ASSAY_USE
  for (j in seq_along(features_list)) { ## this is to reduce num of for loops
    message(features_list[j])
    p3 <- SpatialFeaturePlot(
      vsm_s, features = features_list[j], slot='data', 
      # alpha = c(0.1, 1), 
      pt.size.factor=3) +
      theme(aspect.ratio = get_visium_asp_ratio(vsm_s), legend.position = 'right')  +
      scale_fill_distiller(palette='Spectral', 
                           limits=c(global_feature_value_min[features_list[j]], 
                                    global_feature_value_max[features_list[j]]),
                           oob = scales::squish)
    ggsave2(
      file.path(dir_res, features_list[j], 
                sprintf('SpatialFeaturePlot.%s.%s.%s.%s_%s', 
                        ASSAY_USE, 
                        names(f_vsm_list)[i], features_list[j], 
                        z, paste(z_use, collapse = '_'))),
      p3, width = 5+1, height = 5 * get_visium_asp_ratio(vsm_s))
  }
}
