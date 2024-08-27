#---------------------------
# Compare features across conditions              ----    
#---------------------------
# [ARTC10=ARC2, ARTC23=ARC3]
# Yun Yan
source('util.ruok.R')

dir_res <- './result/viz_features_across_conditions'
fs::dir_create(dir_res)
run_st_list <- c('ART23', 'ART10', 
                 'ART40', 'ART43', ## though CNA looks good, the strength could be a problem
                 'ART65', 'ART304', 
                 'ART305',
                 'ART311', 'ART312')
ASSAY_USE <- 'Spatial'
patient_info <- read_rds('./lib/patient_info.rds')
head(patient_info)
patient_y <- 'pCR'

features_list <- c('HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DMA',
                   'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 
                   'TOP2A', 'MKI67', 'CENPF', 
                   'PCNA', 'TYMS', 'CLSPN', 
                   'ISG15', 'IFI44', 'IFIT3', 
                   'HSPA5', 'DNAJB11', 'SERP1', 
                   'TAGLN', 'ACTA2', 'MYLK', 
                   'S100A6', 'S100A10', 'LGALS3')
features_list <- c('TMSB4X', 'ANXA1', 'KRT19', 'KRT81', 
                   'HIST1H4C', 'PCLAF', 'ATAD2', 'DUT', 'TK1', 'H2AFZ', 
                   'UBE2T', 'TUBB', 'RRM2', 'RAD51AP1', 'DEK')
features_list = c('ANXA1')
f_vsm_list <- file.path('./rds', run_st_list, 'vsm.rds')
names(f_vsm_list) <- run_st_list; all(file.exists(f_vsm_list))

#------ Read ------
## all spots
vsm_list <- lapply(f_vsm_list, read_rds)
## only selected spots
z <- 'confidence'; z_use <- 'High_purity'
print(lapply(vsm_list, function(o) {
  table(o[[z]])
}))
vsm_s_list <- lapply(vsm_list, function(o) {
  Idents(o) <- z
  if (sum(Idents(o) == z_use) == 0) {
    cat('non', z_use, 'spots are available\n')
    # return(subset(o, idents = 'Medium_purity'))
    return(NULL)
  }
  o <- subset(o, idents = z_use)
})
vsm_s_list <- vsm_s_list[!sapply(vsm_s_list, is.null)]
run_st_list <- names(vsm_s_list)
rm(vsm_list)

n_run <- length(run_st_list)
#------ See if gene expression is different across samples ------
# for (j in seq_along(features_list)) {
quick_wilcox <- function(...) {
  o <- wilcox.test(...)
  return(o$p.value)  
}

pal_use <- c('pCR' = '#53D43F', 'non_pCR' = '#811C9A')
pal_pcr <- c('pCR'='#53D43F',
             'RD'='#811C9A', 
             'Unknown'='black', 
             'Excluded'='grey',
             'Removed'='ghostwhite')
for (j in seq_along(features_list)) {
# for (j in c(1)) {  
  g_use <- features_list[j]
  message(g_use)
  fs::dir_create(file.path(dir_res, g_use))
  g_value <- sapply(vsm_s_list, function(o) {
    o@assays[[ASSAY_USE]]@data[g_use, ]
  })
  g_val_df <- enframe_list(g_value, name = 'run', value=g_use)
  pretty_table2str(table(g_val_df$run))
  head(g_val_df)
  g_val_df$y <- tibble::deframe(patient_info[, c('lab_id', patient_y)])[g_val_df$run]
  colnames(g_val_df) <- make.names(colnames(g_val_df))
  p <- ggplot(g_val_df, aes_string(x='run', y=make.names(g_use))) + 
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=2, shape=16, color='peru', alpha=0.5) + 
    stat_compare_means() + 
    scale_x_discrete(labels=pretty_table2str(table(g_val_df$run))) +
    theme_pubr(base_size = 6) +
    theme(axis.ticks.length=unit(0.06, "inch"))
  # print(p)
  # Sys.sleep(0.3)
  ggsave(file.path(dir_res, g_use, sprintf('vlnplot.%s.across_runs.pdf', ASSAY_USE)), 
         p, height = 3, width = 1 * n_run, useDingbats = F)
  
  test_o <- wilcox.test(g_val_df[, make.names(g_use)] ~ g_val_df[, 'y'])
  # boxplot(g_val_df[, make.names(g_use)] ~ g_val_df[, 'y'])
  test_o$p.value
  p <- ggplot(data=g_val_df, aes_string(y='run', x=make.names(g_use))) + 
    
    geom_quasirandom(size=0.8, shape=16, 
                     aes_string(color='y'), 
                     # color='peru', 
                     alpha=0.5) + 
    geom_boxplot(outlier.shape = NA, fill=NA) +
    facet_wrap(~factor(y, levels=c('pCR', 'non_pCR')), nrow=2,
               scales = 'free_y') +
    # facet_wrap(~y, nrow=2, scales = 'free_y') + 
    scale_color_manual(values=pal_use) + 
    labs(caption = sprintf('pval=%s', signif(test_o$p.value, digits = 3))) + 
    theme_pubr(base_size = 6) +
    theme(axis.ticks.length=unit(0.06, "inch")) + 
    rremove('y.title') + rremove('legend')
  # print(p)
  rm(test_o)
  ggsave(file.path(dir_res, g_use, sprintf('vlnplot.%s.across_conditions_%s.test.pdf', ASSAY_USE, patient_y)), 
         p, width = 2, height = 0.6 * n_run, useDingbats = F)
  
 
}

covv <-  c('nCount_Spatial', 'nFeature_Spatial', 'nCount_SCT', 'nFeature_SCT')
for (j in seq_along(covv)) {
  # for (j in c(1)) {  
  covv_use <- covv[j]
  message(covv_use)
  fs::dir_create(file.path(dir_res, covv_use))
  g_value <- sapply(vsm_s_list, function(o) {
    o@meta.data[, covv_use]
  })
  g_val_df <- enframe_list(g_value, name = 'run', value=covv_use)
  pretty_table2str(table(g_val_df$run))
  head(g_val_df)
  g_val_df$y <- tibble::deframe(patient_info[, c('lab_id', patient_y)])[g_val_df$run]
  colnames(g_val_df) <- make.names(colnames(g_val_df))
  p <- ggplot(g_val_df, aes_string(x='run', y=make.names(covv_use))) + 
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=0.8, shape=16, color='peru', alpha=0.5) + 
    stat_compare_means() + 
    scale_x_discrete(labels=pretty_table2str(table(g_val_df$run))) +
    theme_pubr(base_size = 6) +
    theme(axis.ticks.length=unit(0.06, "inch"))
  # print(p)
  # Sys.sleep(0.3)
  ggsave(file.path(dir_res, covv_use, sprintf('vlnplot.%s.across_runs.pdf', covv_use)), 
         p, height = 3, width = 1 * n_run, useDingbats = F)
  
  test_o <- wilcox.test(g_val_df[, make.names(covv_use)] ~ g_val_df[, 'y'])
  p <- ggplot(data=g_val_df, aes_string(y='run', x=make.names(covv_use))) + 
    
    geom_quasirandom(size=0.8, shape=16, 
                     aes_string(color='y'), 
                     # color='peru', 
                     alpha=0.5) + 
    geom_boxplot(outlier.shape = NA, fill=NA) +
    facet_wrap(~factor(y, levels=c('pCR', 'non_pCR')), nrow=2,
               scales = 'free_y') +
    # facet_wrap(~y, nrow=2, scales = 'free_y') + 
    scale_color_manual(values=pal_use) + 
    labs(caption = sprintf('pval=%s', signif(test_o$p.value, digits = 3))) + 
    theme_pubr(base_size = 6) +
    theme(axis.ticks.length=unit(0.06, "inch")) + 
    rremove('y.title') + rremove('legend')
  # print(p)
  rm(test_o)
  ggsave(file.path(dir_res, covv_use, sprintf('vlnplot.%s.across_conditions_%s.test.pdf', covv_use, patient_y)), 
         p, width = 2, height = 0.6 * n_run, useDingbats = F)
  
  
}
