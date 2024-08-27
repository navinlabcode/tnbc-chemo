## Create correlation of all features
##
## Yun Yan
##
source('./ecotype_helper_functions.R')
fs::dir_create('./rfiles')
#------ inputs ------
mat <- read.table('./inputs/Frequency_by_all_cell_states.csv'); dim(mat)
# view(mat)

#------------------ ~~~ start ~~~ --------------------

sample_use <- rownames(mat)
feature_use <- colnames(mat)
igraph_corr_method <- 'spearman'

corr_fea <- cor(mat[sample_use, feature_use, drop=F], 
                method = igraph_corr_method, use = 'pairwise.complete.obs')
corr_fea_pval <- cor_test_pvals_exclude_diag(
  mat[sample_use, feature_use, drop=F], 
  method = igraph_corr_method, 
  p.adjust.method = NULL, #NULL, #'fdr', 
  na.action = 'na.omit')


save.image("./rfiles/ecotype.RData")
