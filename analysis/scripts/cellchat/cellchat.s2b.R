#---------------------------
# bubble plot              ----    
#---------------------------
library(CellChat)
library(tictoc)
library(future)
library(tidyverse)
library(scales)
library(magrittr)
library(ggpubr)
library(patchwork)
dir_res <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/ecotype/cellchat'
cellchat <- read_rds(file.path(dir_res, 'cellchat.rds'))



#------------------- ~~~ pre ~~~ -------------------  
#------ ecotrait membership ------
load(file.path(
  '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102',
  'ecotype', 'use_tumor_hybrid_frac', 'allpatient_with_archetypes', 'use.RData'))
head(memb_feature)
## manually change and match to the main figure
dict_ecotrait_old2new <- structure(
  c(8, 1, 3, 2, 7, 4, 6, 5),
  names = 1:8)
memb_feature <- structure( 
  ruok::replace_vector(memb_feature, dict_ecotrait_old2new),
  names = names(memb_feature) )
memb_feature <- structure(
  paste0('ecotrait', memb_feature),
  names = names(memb_feature) )
## update the feature names
dict_cellstate_old2new <- deframe(unique(cellchat@meta[, c('cellstate_old', 'cellstate')]))
names(memb_feature) <- as.character(dict_cellstate_old2new[names(memb_feature)])
#------ attach 'ecotrait' to cellchat ------
cellchat@meta$ecotrait <- memb_feature[as.character(cellchat@meta$cellstate)]
# ruok::qtable_scatter_hclust(cellchat@meta$ecotrait, cellchat@meta$cellstate)

#------ colors ------
pal_ecotrait <- structure(
  c('#5EB2E3', '#EDE243', '#059C74', '#D26328', 
    '#969595', '#CB78A5', '#0675B2', '#E6A025'), 
  names = paste0('ecotrait', 1:8) ) # matched the main figure
pal_ecotrait_for_cellstate <- structure(
  pal_ecotrait[memb_feature], 
  names = names(memb_feature)
); pal_ecotrait_for_cellstate <- pal_ecotrait_for_cellstate[levels(cellchat@idents)]

#------ util functions ------


#------------------- ~~~ Data trimming ~~~ -------------------  

#------ remove autocrine signaling ------
# ref: https://github.com/sqjin/CellChat/issues/209#issuecomment-851874094
n_cellstates <- dim(cellchat@net$prob)[1]; message(n_cellstates, ' cell states.')
for ( i in 1:n_cellstates ) {
  # print(max(cellchat@net$prob[i, i, ]))
  cellchat@net$prob[i, i, ]   <- 0
  cellchat@net$weight[i, i] <- 0
  cellchat@net$count[i, i]  <- 0
  # cellchat@net$pval[i, i, ]   <- 1
  # print(max(cellchat@net$prob[i, i, ]))
}
#------ remove LR present <= 3 cells ------
cellchat <- filterCommunication(cellchat, min.cells = 3)
df_net_all <- subsetCommunication(cellchat, thresh = 0.05)
LRnum_per_state <- structure(as.numeric(table(cellchat@idents)), names = levels(cellchat@idents))
mat_wt <- cellchat@net$weight
mat_cnt <- cellchat@net$count
# LRnum_per_state <- colMeans(mat_cnt) + rowMeans(mat_cnt)
# identical(names(table(df_net_all$source)), names(table(df_net_all$target)))
# identical(names(table(df_net_all$source)), levels(cellchat@idents))
# LRnum_per_state <- structure(
#   as.numeric(table(df_net_all$source) + table(df_net_all$target)), 
#   names = levels(cellchat@idents))
# barplot(scales::rescale(LRnum_per_state))
# barplot(scales::rescale((LRnum_per_state0)))
# LRnum_per_state <- scales::rescale(LRnum_per_state, to=c(0.1, 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

combo_contri_eco <- read_rds(
  file = file.path(dir_res, 'LR_contribution_groupby_ecotrait.list.rds'))
#------ attach pathway information ------
dict_LR2pathway <- cellchat@LR$LRsig[, c('interaction_name', 'pathway_name')] %>%
  unique() %>% deframe()
df_combo_contri_eco <- do.call('rbind', combo_contri_eco)
all(df_combo_contri_eco$LR %in% rownames(cellchat@LR$LRsig))
df_combo_contri_eco$interaction_name <- df_combo_contri_eco$LR
df_combo_contri_eco$pathway_name <- cellchat@LR$LRsig[df_combo_contri_eco$LR, 'pathway_name']
df_combo_contri_eco$interaction_name_2 <- df_combo_contri_eco$name

colnames(df_combo_contri_eco)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#------ bubble plot for ecotrait without considering ecotrait specificity ------
dir_viz <- file.path(dir_res, 'netVisual_bubble_ecotrait_global')
fs::dir_create(dir_viz)
ecotrait_focus = 'ecotrait1'
for (ecotrait_focus in names(pal_ecotrait)) {
  message(ecotrait_focus)
  cellstates_focus <- dplyr::filter(cellchat@meta, ecotrait %in% ecotrait_focus) %>%
    dplyr::pull(cellstate) %>% unique()
  print(cellstates_focus)
  
  p <- netVisual_bubble(
    cellchat, 
    sources.use = cellstates_focus,
    targets.use = cellstates_focus,
    remove.isolate = FALSE, 
    return.data = T)
  class(p)
  pmat <- p$communication
  colnames(pmat)
  pmat <- pivot_wider(
    pmat[, c('source.target', 'interaction_name_2', 'prob')], 
    names_from = interaction_name_2, 
    values_from = prob)
  dict_src_tar_labels <- structure(
    pmat$source.target, 
    names = make.names(pmat$source.target)
  )
  pmat <- pmat %>% dplyr::select(-source.target)
  pmat <- as.matrix(pmat)
  rownames(pmat) <- names(dict_src_tar_labels)
  dim(pmat)
  pmat[1:3, 1:3]
  pmat[is.na(pmat)] <- 0
  hclust_ST <- hclust(dist(pmat))
  hclust_LR  <- hclust(dist(t(pmat)))
  
  
  pnew <- p$gg.obj + 
    scale_x_discrete(limits = as.character(dict_src_tar_labels[hclust_ST$order])) + 
    scale_y_discrete(limits = as.character(colnames(pmat)[hclust_LR$order]) )
  
  pdf(file.path(dir_viz, sprintf('%s.bubble.pdf', ecotrait_focus)), 
      width = 20, height = 30)
  print(pnew)
  dev.off()
}
  
#------ bubble plot for each ecotrait with only cancer~TME communications ------
dict_cellstate_is_cancer <- structure(
  rep(FALSE, length(names(memb_feature))), 
  names = names(memb_feature)
)
dict_cellstate_is_cancer[1:11] <- T

for (ecotrait_focus in names(pal_ecotrait)) {
  message(ecotrait_focus)
  cellstates_focus <- dplyr::filter(cellchat@meta, ecotrait %in% ecotrait_focus) %>%
    dplyr::pull(cellstate) %>% unique()
  print(cellstates_focus)
  
  p <- netVisual_bubble(
    cellchat, 
    sources.use = cellstates_focus,
    targets.use = cellstates_focus,
    remove.isolate = FALSE, 
    return.data = T)
  class(p)
  pmat <- p$communication
  colnames(pmat)
  head(pmat)
  all(pmat$target %in% names(dict_cellstate_is_cancer))
  all(pmat$source %in% names(dict_cellstate_is_cancer))
  
  i1 <- dict_cellstate_is_cancer[as.character(pmat$source)] & !dict_cellstate_is_cancer[as.character(pmat$target)] ; sum(i1)
  i2 <- !dict_cellstate_is_cancer[as.character(pmat$source)] & dict_cellstate_is_cancer[as.character(pmat$target)]; sum(i2)
  sum(i1 | i2)
  pmat <- pmat[which(i1|i2), , drop=F]
  if (nrow(pmat) == 0) {warning('no cancer - TME LR found'); next()}
  
  pmat <- pivot_wider(
    pmat[, c('source.target', 'interaction_name_2', 'prob')], 
    names_from = interaction_name_2, 
    values_from = prob)
  dict_src_tar_labels <- structure(
    pmat$source.target, 
    names = make.names(pmat$source.target)
  )
  pmat <- pmat %>% dplyr::select(-source.target)
  pmat <- as.matrix(pmat)
  rownames(pmat) <- names(dict_src_tar_labels)
  dim(pmat)
  # pmat[1:3, 1:3]
  pmat[is.na(pmat)] <- 0
  hclust_ST <- hclust(dist(pmat))
  hclust_LR  <- hclust(dist(t(pmat)))
  
  n_ST <- nrow(pmat)
  n_LR <- ncol(pmat)
  pnew <- p$gg.obj + 
    scale_x_discrete(limits = as.character(dict_src_tar_labels[hclust_ST$order])) + 
    scale_y_discrete(limits = as.character(colnames(pmat)[hclust_LR$order]) ) 
  
  pdf(file.path(dir_viz, sprintf('%s.bubble.cancer_TME.pdf', ecotrait_focus)), 
      width = pmax(n_ST * 0.2, 5), height = n_LR * 0.15)
  print(pnew)
  dev.off()
}
  
  