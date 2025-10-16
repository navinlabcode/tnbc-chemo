#---------------------------
# Find LR that are intra-ecotrait-specific:
# high probs inside an ecotrait and low probs outside ecotrait              ----    
# 
# But quite hard. 
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
write_rds(cellchat, file.path(dir_res, 'cellchat.original.rds'))

#------------------- ~~~ pre ~~~ -------------------  
#------ ecotrait membership ------
load(file.path('/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102',
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

#------------------- ~~~ Viz global interactions ~~~ -------------------  
pdf(file.path(dir_res, 'netVisual_heatmap.global_count.pdf'), width = 8, height = 7, onefile = T)
netVisual_heatmap(
  cellchat, measure = 'count', 
  color.use = pal_ecotrait_for_cellstate,
  color.heatmap = "Reds")
netVisual_heatmap(
  cellchat, measure = 'weight', 
  color.use = pal_ecotrait_for_cellstate,
  color.heatmap = "Purples")
dev.off()

#------------------- ~~~ Viz weight with each ecotrait ~~~ -------------------  
mat_wt <- cellchat@net$weight
ecotrait_opts <- names(pal_ecotrait)
mat_wt[mat_wt < quantile(mat_wt, 0.5)] <- 0
pdf(file.path(dir_res, 'netVisual_circle.weight_per_ecotrait.pdf'), width = 7, height = 7, onefile = T, useDingbats = F)
for (ecotrait_focus in ecotrait_opts) {
  message(ecotrait_focus)
  mat2 <- matrix(0, nrow = nrow(mat_wt), ncol = ncol(mat_wt), 
                 dimnames = dimnames(mat_wt))
  cellstates_focus <- names(memb_feature)[memb_feature == ecotrait_focus]
  mat2[cellstates_focus, ] <- mat_wt[cellstates_focus, ]
  mat2[, cellstates_focus] <- mat_wt[, cellstates_focus]
  netVisual_circle(
    mat2,
    edge.weight.max = max(mat_wt), 
    title.name = ecotrait_focus, 
    vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F,
    edge.width.max = 1,
    color.use = pal_ecotrait_for_cellstate)
  
}
dev.off()
#------------------- ~~~ Evaluate LR contribution for intra-/inter- ecotrait ~~~ -------------------  
#------ calc LR contributions in intra- and inter-ecotrait ------
# // ST = LR contribution intra-ecotrait
# // S  = LR contribution of inter-ecotrait as sender
# // T  = LR contribution of inter-ecotrait as target
# // intra-ecotrait contribution = ST
# // inter-ecotrait contribution = ((S - ST) + (T - ST)) / 2
# 
# LR ~ contri_intra_eco
# LR ~ contri_inter_eco
# Prefer LR which has contri_intra_eco > contri_inter_eco. 

## viz
combo_contri_eco <- lapply(names(pal_ecotrait), function(ecotrait_focus) {
  
  message(ecotrait_focus)
  
  cellstates_focus <- dplyr::filter(cellchat@meta, ecotrait %in% ecotrait_focus) %>%
    dplyr::pull(cellstate) %>% unique()
  print(cellstates_focus)
  
  if (FALSE) {
    netVisual_circle(
      cellchat@net$weight,

      vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F,
      color.use = pal_ecotrait_for_cellstate)

    netVisual_circle(
      cellchat@net$weight, 
      sources.use = cellstates_focus, 
      
      vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F, 
      color.use = pal_ecotrait_for_cellstate)
    
    netVisual_circle(
      cellchat@net$weight, 
      targets.use = cellstates_focus, 
      
      vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F, 
      color.use = pal_ecotrait_for_cellstate)
    
    netVisual_circle(
      cellchat@net$weight, 
      sources.use = cellstates_focus, targets.use = cellstates_focus, 
      
      vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F, 
      color.use = pal_ecotrait_for_cellstate)
  }
  
  cat('---- extract LR data\n')
  if (F) {
    contri_eco_src <- try(CellChat::netAnalysis_contribution_allLR(
      cellchat, return.data = T, 
      sources.use = as.character(cellstates_focus), 
      targets.use = levels(cellchat@idents)))
    contri_eco_tar <-  try(CellChat::netAnalysis_contribution_allLR(
      cellchat, return.data = T, 
      sources.use = levels(cellchat@idents), 
      targets.use = as.character(cellstates_focus)))
    contri_intra_eco <- try( CellChat::netAnalysis_contribution_allLR(
      cellchat, return.data = T, 
      sources.use = as.character(cellstates_focus), 
      targets.use = as.character(cellstates_focus)))
  } else {
    contri_eco_src <- try(netAnalysis_contribution_allLR(
      cellchat, return.data = T, 
      sources.use = as.character(cellstates_focus), 
      targets.use = levels(cellchat@idents)))
    contri_eco_tar <-  try(netAnalysis_contribution_allLR(
      cellchat, return.data = T, 
      sources.use = levels(cellchat@idents), 
      targets.use = as.character(cellstates_focus)))
    contri_intra_eco <- try(netAnalysis_contribution_allLR(
      cellchat, return.data = T, 
      sources.use = as.character(cellstates_focus), 
      targets.use = as.character(cellstates_focus)))
  }
  
  contri_eco_src   <- contri_eco_src$LR.contribution
  contri_eco_tar   <- contri_eco_tar$LR.contribution
  contri_intra_eco <- contri_intra_eco$LR.contribution
  
  nrow(contri_eco_src)
  nrow(contri_eco_tar)
  nrow(contri_intra_eco)
  
  # df_net_intra_eco <- subsetCommunication(
  #   cellchat, 
  #   targets.use = as.character(cellstates_focus),
  #   sources.use = as.character(cellstates_focus))
  
  contri_eco_src$LR <- rownames(contri_eco_src)
  contri_eco_tar$LR <- rownames(contri_eco_tar)
  contri_intra_eco$LR <- rownames(contri_intra_eco)
  all(contri_intra_eco$LR %in% contri_eco_src$LR)
  all(contri_intra_eco$LR %in% contri_eco_tar$LR)
  
  cat('---- calculate the pure inter-ecotrait interaction.\n')
  contri_inter_src <- dplyr::left_join(
    x=contri_eco_src, y = contri_intra_eco,
    by = c('LR', 'name'), suffix = c('_src_all', '_intra'))
  
  contri_inter_tar <- dplyr::left_join(
    x=contri_eco_tar, y = contri_intra_eco,
    by = c('LR', 'name'), suffix = c('_tar_all', '_intra'))
  
  replace_na <- function(x, na_to=0) { x[is.na(x)] <- na_to; return(x) }
  contri_inter_src$contribution_intra <- replace_na(contri_inter_src$contribution_intra)
  contri_inter_tar$contribution_intra <- replace_na(contri_inter_tar$contribution_intra)
  
  
  contri_inter_src %<>% dplyr::mutate(
    contribution_src = contribution_src_all - contribution_intra)
  contri_inter_tar %<>% dplyr::mutate(
    contribution_tar = contribution_tar_all - contribution_intra)
  
  contri_inter_src %<>% dplyr::select(-contribution_intra)
  contri_inter_tar %<>% dplyr::select(-contribution_intra)
  
  contri_inter_eco <- dplyr::full_join(
    x=contri_inter_src, y=contri_inter_tar, 
    by = c('LR', 'name'), 
    suffix = c('_SRC', '_TAR'))
  colnames(contri_inter_eco)
  contri_inter_eco$contribution_src <- replace_na(contri_inter_eco$contribution_src)
  contri_inter_eco$contribution_tar <- replace_na(contri_inter_eco$contribution_tar)
  contri_inter_eco %<>% dplyr::mutate(
    contribution = (contribution_src + contribution_tar)/2)
  contri_inter_eco %<>% dplyr::select(
    -contribution_src_all, -contribution_tar_all
  )
  
  ## compare to find LR specific to eco
  cat('---- rank the LR with intra>inter.\n')
  head(contri_intra_eco)
  head(contri_inter_eco)
  
  contri_eco <- dplyr::full_join(
    x=contri_intra_eco, y=contri_inter_eco, 
    by = c('LR', 'name'), 
    suffix = c('_intra', '_inter')
  )
  
  contri_eco %<>% dplyr::mutate(
    contribution_intra = tidyr::replace_na(contribution_intra, 0), 
    contribution_src = tidyr::replace_na(contribution_src, 0), 
    contribution_tar = tidyr::replace_na(contribution_tar, 0), 
    contribution_inter = tidyr::replace_na(contribution_inter, 0)
  )
  contri_eco %<>% dplyr::mutate(
    diff_intra_inter = contribution_intra - contribution_inter
  )
  contri_eco$ecotrait <- ecotrait_focus
  contri_eco
  ## finish each ecotrait
})
names(combo_contri_eco) <- names(pal_ecotrait)
print(sapply(combo_contri_eco, nrow))


write_rds(
  combo_contri_eco, 
  file = file.path(dir_res, 'LR_contribution_groupby_ecotrait.list.rds'))
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


#------------------- ~~~ Propose LR simply by ranking top [deprecated] ~~~ -------------------  
# ecotrait1:
#   path1, path2
#   {lr1, lr2}, {lr3, lr4}
# ecotrait2:
#   path1, path3
#   {lr10, lr20}, {lr5, lr6}
# pathway-centric analysis would have problem. 
# todo!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
quantile(df_combo_contri_eco$diff_intra_inter)
quantile(df_combo_contri_eco$contribution_intra)
colnames(df_combo_contri_eco)
mat_contri_LR_eco <- df_combo_contri_eco %>%
  dplyr::filter(contribution_intra >= quantile(contribution_intra, 0.9)) %>%
  # dplyr::filter(diff_intra_inter>0) %>%
  # dplyr::filter(diff_intra_inter>=quantile(df_combo_contri_eco$diff_intra_inter,.1)) %>%
  # dplyr::filter(diff_intra_inter>0, 
  #               contribution_intra>contribution_src,contribution_intra>contribution_tar) %>%
  # dplyr::select(interaction_name, ecotrait, diff_intra_inter) %>%
  # tidyr::pivot_wider(names_from = ecotrait, values_from = diff_intra_inter)
  dplyr::select(interaction_name, ecotrait, contribution_intra) %>%
  tidyr::pivot_wider(names_from = ecotrait, values_from = contribution_intra)


mat_contri_LR_eco <- as.data.frame(mat_contri_LR_eco)
rownames(mat_contri_LR_eco) <- mat_contri_LR_eco$interaction_name
mat_contri_LR_eco$interaction_name <- NULL
mat_contri_LR_eco <- as.matrix(mat_contri_LR_eco)
mat_contri_LR_eco[is.na(mat_contri_LR_eco)] <- -Inf
dict_highestLR_eco <- apply(mat_contri_LR_eco, 1, nnet::which.is.max)
dict_highestLR_eco <- structure(
  colnames(mat_contri_LR_eco)[dict_highestLR_eco], 
  names = names(dict_highestLR_eco)
)
table(dict_highestLR_eco)
# ecotrait1 ecotrait2 ecotrait3 ecotrait4 
# 35        41        61       115 
# ecotrait5 ecotrait6 ecotrait7 ecotrait8 
# 60        13        71        53
# ecotrait3 ecotrait4 ecotrait5 ecotrait6 
# 5         5         3         1 
# ecotrait7 ecotrait8 
# 1         2 

## for ecotrait1, find some LR from the 41 LRs for visualization
df_highestLR_eco <- enframe(
  dict_highestLR_eco, name = 'interaction_name', value = 'ecotrait') %>%
  as.data.frame()
rownames(df_highestLR_eco) <- df_highestLR_eco$interaction_name
head(df_highestLR_eco)
df_highestLR_eco$strength <- sapply(1:nrow(df_highestLR_eco), function(r) {
  i <- df_highestLR_eco[r, 'interaction_name']
  j <- df_highestLR_eco[r, 'ecotrait']
  mat_contri_LR_eco[i, j]
})
df_highestLR_eco %<>% dplyr::arrange(ecotrait, desc(strength))
table(df_highestLR_eco$ecotrait)

view(df_highestLR_eco)

LR_i <- 'CD22_PTPRC'
fs::dir_create(file.path(dir_res, 'netVisual_individual_LR'))
for (LR_i in seq_len(nrow(df_highestLR_eco))) {
  cat(LR_i, ' ')
  ecotrait_focus <- df_highestLR_eco[LR_i, 'ecotrait']
  message(ecotrait_focus)
  cellstates_focus <- dplyr::filter(
    cellchat@meta, ecotrait %in% ecotrait_focus) %>%
    dplyr::pull(cellstate) %>% unique()
  
  # mat_LR_i <- cellchat@net$prob[, , df_highestLR_eco[LR_i, 'interaction_name']]
  # library(ComplexHeatmap)
  # Heatmap(as.matrix(mat_LR_i>0)*1)
  # netVisual_circle(mat_LR_i)
  pdf(file.path(dir_res, 'netVisual_individual_LR', 
                sprintf('%s.LR_%s.pdf', ecotrait_focus, 
                        df_highestLR_eco[LR_i, 'interaction_name'])), 
      width = 7, height = 7, useDingbats = F)
  netVisual_individual(
    cellchat,
    color.use = pal_ecotrait_for_cellstate,
    # group = memb_feature[levels(cellchat@idents)],
    signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
    pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
    remove.isolate = F,
    # vertex.size.max = 0.01,
    # vertex.weight = 1, 
    layout = 'chord',
    # vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F
  ) 
  dev.off()
}  

netAnalysis_contribution_allLR(
  cellchat, 
  sources.use = as.character(cellstates_focus), 
  targets.use = as.character(cellstates_focus))

netVisual_bubble(
  cellchat, 
  sources.use = cellstates_focus, 
  targets.use = cellstates_focus, remove.isolate = T)


netVisual_chord_gene(
  cellchat, 
  # sources.use = cellstates_focus, targets.use = cellstates_focus,
  pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F]
)

netVisual_chord_gene(
  cellchat, 
  sources.use = cellstates_focus, targets.use = cellstates_focus,
  pairLR.use = df_highestLR_eco[16:17, 'interaction_name', drop=F]
)
netVisual_chord_gene(
  cellchat, 
  sources.use = cellstates_focus, targets.use = cellstates_focus,
  # pairLR.use = df_highestLR_eco[16:17, 'interaction_name', drop=F]
)


netVisual_individual(
  cellchat,
  sources.use = cellstates_focus, 
  targets.use = cellstates_focus,
  color.use = pal_ecotrait_for_cellstate,
  signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
  pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
  remove.isolate = F,
  vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F
)
netVisual_individual(
  cellchat,
  sources.use = cellstates_focus, 
  
  color.use = pal_ecotrait_for_cellstate,
  signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
  pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
  remove.isolate = F,
  vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F
)
netVisual_individual(
  cellchat,

  targets.use = cellstates_focus,
  color.use = pal_ecotrait_for_cellstate,
  signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
  pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
  remove.isolate = F,
  vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F
)


#------------------- ~~~ Propose ecotrait-specific LR ~~~ -------------------  
# Using the intra_ecotrait contribution, find LR that is highly specific to ecotrait
# caveat: does not ensure the LR has no/smaller inter-ecotrait communications

dict_interaction_name2bio <- cellchat@LR$LRsig[, c('interaction_name', 'interaction_name_2')] %>% deframe()
#------ consider the intra-eco contribution ------
mat_contri_LR_eco <- df_combo_contri_eco %>%
  # dplyr::filter(contribution_intra >= quantile(contribution_intra, 0.9)) %>%
  # dplyr::filter(diff_intra_inter>0) %>%
  # dplyr::filter(diff_intra_inter>=quantile(df_combo_contri_eco$diff_intra_inter,.1)) %>%
  # dplyr::filter(diff_intra_inter>0, 
  #               contribution_intra>contribution_src,contribution_intra>contribution_tar) %>%
  # dplyr::select(interaction_name, ecotrait, diff_intra_inter) %>%
  # tidyr::pivot_wider(names_from = ecotrait, values_from = diff_intra_inter)
  dplyr::select(interaction_name, ecotrait, contribution_intra) %>%
  tidyr::pivot_wider(names_from = ecotrait, values_from = contribution_intra)
mat_contri_LR_eco <- as.data.frame(mat_contri_LR_eco)
rownames(mat_contri_LR_eco) <- mat_contri_LR_eco$interaction_name
mat_contri_LR_eco$interaction_name <- NULL
mat_contri_LR_eco <- as.matrix(mat_contri_LR_eco)
mat_contri_LR_eco[is.na(mat_contri_LR_eco)] <- -Inf

#------ find the marker LR per ecotrait ------
  
## apply the same strategy to find the marker LR per ecotrait
## re-use the function
define_marker_gene_per_nsnmf <- function(X, n_tolerance=0){
  # X: NMF gene loading matrix (W) genes x factors
  # Ref: Reuben and Itai Yanai's ST paper.
  ## Don't use rank=1 because of ties. 
  factor_names <- colnames(X)
  mat_row_ismax <- t( apply(X, 1, function(xx) xx == max(xx, na.rm=TRUE)) )
  colnames(mat_row_ismax) <- factor_names
  res <- lapply(factor_names, function(J){
    n_violate <- 0
    v <- sort(X[, J], decreasing = TRUE, na.last=TRUE)
    gnames <- names(v)
    ismax <- mat_row_ismax[gnames, J]
    out <- c()
    i <- 1
    while( i <= length(gnames) ){
      if (!ismax[i]) { n_violate <- n_violate + 1}
      if (n_violate > n_tolerance) { break() }
      if (ismax[i] ) {out <- c(out, gnames[i])}
      i <- i + 1
    }
    return(out)
  })
  names(res) <- factor_names
  return(res)
}
nrow(mat_contri_LR_eco)

list_topLR_eco <- define_marker_gene_per_nsnmf(
  mat_contri_LR_eco, n_tolerance = 100)
list_topLR_eco <- lapply(list_topLR_eco, head, 10)
list_topLR_eco
table(duplicated(unlist(list_topLR_eco)))
df_highestLR_eco <- ruok::enframe_list(list_topLR_eco, 'ecotrait', 'interaction_name')
rownames(df_highestLR_eco) <- df_highestLR_eco$interaction_name
table(duplicated(df_highestLR_eco$interaction_name))

df_highestLR_eco$strength <- sapply(1:nrow(df_highestLR_eco), function(r) {
  i <- df_highestLR_eco[r, 'interaction_name']
  j <- df_highestLR_eco[r, 'ecotrait']
  mat_contri_LR_eco[i, j]
})
df_highestLR_eco %<>% dplyr::arrange(ecotrait, desc(strength))
df_highestLR_eco$plot_order <- nrow(df_highestLR_eco):1
head(df_highestLR_eco)

df_highestLR_eco$interaction_name_bio <- dict_interaction_name2bio[df_highestLR_eco$interaction_name]
#------ absolute prob intra-ecotrait ------
pa <- ggplot(df_highestLR_eco, 
             aes(x=strength, y=interaction_name)) + 
  geom_col(aes(fill=ecotrait)) +
  scale_fill_manual(values = pal_ecotrait) +
  scale_y_discrete(limits = df_highestLR_eco$interaction_name[df_highestLR_eco$plot_order]) +
  theme_bw()
pa
pa <- ggplot(df_highestLR_eco, 
             aes(x=strength, y=reorder(interaction_name_bio, plot_order))) + 
  geom_col(aes(fill=ecotrait)) +
  scale_fill_manual(values = pal_ecotrait) +
  facet_wrap(~ecotrait, ncol = 1, scales = 'free', strip.position='left') + 
  theme_bw(base_size = 4)
pa <- pa + rremove('legend')
pa
#------ relative prob inter-ecotrait per LR ------
df_rel_prop_per_LR <- df_combo_contri_eco %>%
  dplyr::filter(interaction_name %in% unique(df_highestLR_eco$interaction_name)) %>%
  dplyr::mutate(strength = contribution_intra) %>%
  dplyr::mutate(topecotrait = deframe(df_highestLR_eco[, c('interaction_name', 'ecotrait')])[.data[['interaction_name']]]) %>%
  dplyr::group_by(interaction_name) %>%
  dplyr::mutate(strength = rescale(strength)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(interaction_name_bio = dict_interaction_name2bio[.data[['interaction_name']]]) %>%
  dplyr::mutate(fill_col = ifelse(ecotrait == topecotrait, .data[['topecotrait']], 'other'))
pb <- df_rel_prop_per_LR %>%
  ggplot(aes(x=strength, 
             y=interaction_name_bio, 
             fill=fill_col)) + 
  geom_col() + 
  facet_wrap(~ecotrait, nrow = 1) +
  scale_y_discrete(
    limits = df_highestLR_eco$interaction_name_bio[df_highestLR_eco$plot_order]) +
  # scale_fill_manual(values = c(`TRUE`='orange', `FALSE`='black')) + 
  scale_fill_manual(values = c(pal_ecotrait, 'other'='black')) + 
  labs(x='relative strength') +
  theme_bw(base_size = 4) + rremove('legend')


pdf(file.path(dir_res, 'barplot.markerLR_per_ecotrait.pdf'), 
    width = 7, height = 5, useDingbats = F)
print(wrap_plots(pa, pb, widths = c(1, 3)))
dev.off()

pb_heatmap <- df_rel_prop_per_LR %>%
  dplyr::select(interaction_name, ecotrait, strength) %>%
  tidyr::pivot_wider(names_from = ecotrait, values_from = strength) %>%
  as.data.frame()
rownames(pb_heatmap) <- pb_heatmap$interaction_name
pb_heatmap$interaction_name <- NULL
pb_heatmap <- as.matrix(pb_heatmap)
pb_heatmap[is.na(pb_heatmap)] <- 0
pb_heatmap <- pb_heatmap[df_highestLR_eco$interaction_name, ]
library(ComplexHeatmap)
ht <- Heatmap(
  pb_heatmap, name = 'relative strength',
  col = circlize::colorRamp2(c(0, 0.5, 1), c('ghostwhite', 'yellow', 'red')),
  cluster_rows = F, cluster_columns = F,
  show_row_names = T, show_column_names = T,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = 'relative strength',
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 6)
  )
)
pdf(file.path(dir_res, 'heatmap.markerLR_per_ecotrait.pdf'), 
    width = 5, height = 7, useDingbats = F)
draw(ht)
dev.off()

#------ viz each LR in graph ------
list_topLR_eco
fs::dir_create(file.path(dir_res, 'netVisual_individual_LR_new'))
for (LR_i in 1:nrow(df_highestLR_eco)) {
# for (LR_i in c(1)) {
  message(LR_i, '/', nrow(df_highestLR_eco))
  ecotrait_focus <- df_highestLR_eco[LR_i, 'ecotrait']
  message(ecotrait_focus, ' ', df_highestLR_eco[LR_i, 'interaction_name'])
  cellstates_focus <- dplyr::filter(
    cellchat@meta, ecotrait %in% ecotrait_focus) %>%
    dplyr::pull(cellstate) %>% unique()

  pdf(file.path(dir_res, 'netVisual_individual_LR_new',
                paste0(sprintf('%s.LR_%s', ecotrait_focus,
                               make.names(df_highestLR_eco[LR_i, 'interaction_name'])), 
                       '%2d.pdf')
                ),
      width = 5, height = 5, useDingbats = F, onefile = F)
  # all cell states
  netVisual_individual(
    cellchat,
    color.use = pal_ecotrait_for_cellstate,
    group = memb_feature[levels(cellchat@idents)],
    # sources.use = cellstates_focus, targets.use = cellstates_focus,
    signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
    pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
    remove.isolate = T,
    layout = 'chord')
  # intra-eco cell states
  netVisual_individual(
    cellchat,
    # color.use = pal_ecotrait_for_cellstate,
    # group = memb_feature[levels(cellchat@idents)],
    sources.use = cellstates_focus, targets.use = cellstates_focus,
    signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
    pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
    remove.isolate = T,
    layout = 'circle')
  dev.off()
}
fs::dir_create(file.path(dir_res, 'netVisual_bubble_ecotrait'))
for (ecotrait_focus in unique(df_highestLR_eco$ecotrait) ) {
  # for (LR_i in c(1)) {
  # message(LR_i, '/', nrow(df_highestLR_eco))
  # ecotrait_focus <- df_highestLR_eco[LR_i, 'ecotrait']
  # message(ecotrait_focus, ' ', df_highestLR_eco[LR_i, 'interaction_name'])
  
  
  cellstates_focus <- dplyr::filter(
    cellchat@meta, ecotrait %in% ecotrait_focus) %>%
    dplyr::pull(cellstate) %>% unique()
  pdf(file.path(dir_res, 'netVisual_bubble_ecotrait',
                sprintf('netVisual_bubble_ecotrait.%s.pdf', ecotrait_focus)),
  width = 10, height = 10, useDingbats = F, onefile = F)
  
  p <- netVisual_bubble(
    cellchat, 
    sources.use = cellstates_focus,
    targets.use = cellstates_focus,
    # signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
    pairLR.use = df_highestLR_eco[df_highestLR_eco$ecotrait %in% ecotrait_focus, 'interaction_name', drop=F], 
    remove.isolate = FALSE, 
    return.data = F)
  print(p)
  dev.off()
  

}


# netVisual_chord_gene(
#   cellchat, 
#   sources.use = cellstates_focus, targets.use = cellstates_focus,
#   pairLR.use = df_highestLR_eco[list_topLR_eco$ecotrait8, 'interaction_name', drop=F]
# )
# netVisual_chord_gene(
#   cellchat, 
#   sources.use = cellstates_focus, targets.use = cellstates_focus,
#   pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F]
# )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------- ~~~ Supervised way to zoom in certain cell states ~~~ -------------------
pal_ecotrait_for_cellstate
cellstates_src <- c('M05__Interferon')
cellstates_tar <- c('CD4-TIFN', 'BIFN', 'CD8-TIFN', 'mac-IFN')
netVisual_individual
netVisual_chord_gene(
  cellchat,
  # color.use = pal_ecotrait_for_cellstate,
  # group = memb_feature[levels(cellchat@idents)],
  # sources.use = cellstates_src, targets.use = cellstates_tar,
  # sources.use = cellstates_tar, targets.use = cellstates_src,
  signaling = c('GALECTIN', 'COMPLEMENT'),
  # signaling = dict_LR2pathway[df_highestLR_eco[LR_i, 'interaction_name']],
  # pairLR.use = df_highestLR_eco[LR_i, 'interaction_name', drop=F], 
  # remove.isolate = F,
  # top=0.5,
  # slot.name = "netP",
  # layout = 'circle'
  )

#------------------- ~~~ [deprecated below] ~~~ -------------------  


#------ report top 3 LRs per pathway per ecotrait ------
list_combo_contri_eco <- split(df_combo_contri_eco, df_combo_contri_eco$ecotrait)
names(list_combo_contri_eco)

find_unique_per_slot <- function(l){
  o <- sapply(1:length(l), function(i) {
    v_i <- l[[i]]
    v_j <- unlist(l[setdiff(1:length(l), i)])
    return(setdiff(v_i, v_j))
  })
  names(o) <- names(l)
  return(o)
}

ecotrait_specific_pathways <- sapply(
  list_combo_contri_eco, 
  function(df) {
    dplyr::filter(df, diff_intra_inter > 0) %>%
      dplyr::pull(pathway_name) %>%
      unique()
  }
)
sapply(ecotrait_specific_pathways, length)

ecotrait_exclusive_pathways <- find_unique_per_slot(ecotrait_specific_pathways)
sapply(ecotrait_exclusive_pathways, length)

ecotrait_focus <- 'ecotrait6'
cellstates_focus <- dplyr::filter(
  cellchat@meta, ecotrait %in% ecotrait_focus) %>%
  dplyr::pull(cellstate) %>% unique()

df_combo_contri_eco_viz <- df_combo_contri_eco %>% 
  dplyr::filter(diff_intra_inter > 0, ecotrait %in% ecotrait_focus) %>%
  dplyr::group_by(ecotrait, pathway_name) %>% 
  # dplyr::slice_max(diff_intra_inter, n=3, with_ties = F) %>%
  dplyr::arrange(desc(diff_intra_inter))

netVisual_circle(
  cellchat@net$weight, 
  sources.use = cellstates_focus, 
  targets.use = cellstates_focus,
  vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F, 
  color.use = pal_ecotrait_for_cellstate)

df_combo_contri_eco_viz <- df_combo_contri_eco %>% 
  dplyr::filter(pathway_name %in% ecotrait_exclusive_pathways[[ecotrait_focus]]) %>%
  dplyr::arrange(desc(diff_intra_inter))
netVisual_individual(
  cellchat,
  # sources.use = cellstates_focus, 
  # targets.use = cellstates_focus,
  # color.use = pal_ecotrait_for_cellstate[cellstates_focus], 
  color.use = pal_ecotrait_for_cellstate,
  signaling = df_combo_contri_eco_viz[1, 'pathway_name', drop=T],
  pairLR.use = df_combo_contri_eco_viz[1, 'interaction_name', drop=F], 
  remove.isolate = F, 
  edge.curved = 0.5,
  vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F
)

netVisual_individual(
  cellchat,
  sources.use = cellstates_focus,
  targets.use = cellstates_focus,
  color.use = pal_ecotrait_for_cellstate,
  signaling = df_combo_contri_eco_viz[1, 'pathway_name', drop=T],
  pairLR.use = df_combo_contri_eco_viz[1, 'interaction_name', drop=F], 
  # vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F, 
  remove.isolate = F
)
netVisual_circle(
  cellchat@net$weight, 
  sources.use = cellstates_focus, targets.use = cellstates_focus, 
  
  vertex.weight = LRnum_per_state, weight.scale = T, label.edge= F, 
  color.use = pal_ecotrait_for_cellstate)

netVisual_heatmap(
  cellchat, 
  # sources.use = cellstates_focus,
  # targets.use = cellstates_focus,
  signaling = df_combo_contri_eco_viz[1, 'pathway_name', drop=T],
  color.use = pal_ecotrait_for_cellstate,
  color.heatmap = "Reds")

df_net_intra_eco <- subsetCommunication(
  cellchat, 
  sources.use = cellstates_focus, 
  targets.use = cellstates_focus)

