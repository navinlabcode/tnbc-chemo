#---------------------------
# Graph similarity of pCR/RD and Archetypes
# Metric: jaccard similarity             ----   
## Yun Yan 
#---------------------------

if (T) {
  library(matrixStats)
  library(Matrix)
  library(ggpubr); library(igraph); library(tidyverse)
  library(scales)
  source('util.nmf.viz.R')
}

dir_proj <- file.path(
  '.',
  'ecotype', 'use_tumor_hybrid_frac')


load("./rfiles/ecotype.RData")


dir_res <- file.path(dir_proj, 'feature_galaxy.copresence.by_obs_groups')
fs::dir_create(dir_res)

value_type <- 'zscore'
z_type_opts <- c('pCR_status', 'archetype')

#------ read the original raw graph ------
  
graph_combo_response <- combo_response

graph_combo_archetype <- combo_archetypes


graph_combo_response <- lapply(graph_combo_response, function(obj) {
  Gf <- obj$graph
  bad_edge <- (! E(Gf)$label %in% igraph_feature_vizsig_edges)
  bad_edge[is.na(bad_edge)] <- T
  return( delete_edges(Gf, E(Gf)[bad_edge]) )  
})
graph_combo_archetype <- lapply(graph_combo_archetype, function(obj) {
  Gf <- obj$graph
  bad_edge <- (! E(Gf)$label %in% igraph_feature_vizsig_edges)
  bad_edge[is.na(bad_edge)] <- T
  return( delete_edges(Gf, E(Gf)[bad_edge]) )  
})



## function for jaccard similarity of edge sets
jaccard_edgeset_similarity <- function(G1, G2) {
  inter <- length(E(G1 %s% G2))
  un <- length(E(G1 %u% G2))
  
  if (un == 0) {
    0
  } else {
    inter/un
  }
}


### Reponse vs Archetypes
mat_graph_simi <- sapply(graph_combo_archetype, function(a) {
  sapply(graph_combo_response, function(b){
    jaccard_edgeset_similarity(a, b)
  })
})

set_aspect_ratio <- function(mat, ratio=1) {
  ncol(mat) / nrow(mat) * ratio
}

pdf(file.path(dir_res, sprintf('heatmap_graph_similarity_response_x_archetype_%s.pdf', value_type)), 
    width = 1.5*ncol(mat_graph_simi)/nrow(mat_graph_simi) + 1, 
    height = 1.5+1, 
    useDingbats = F)
p <- Heatmap(mat_graph_simi, name = 'jaccard', 
             row_names_side = 'left', column_names_side = 'top',
             cluster_rows = F, cluster_columns = F, 
             col = heatmap_color_fun_cont(from=0, 
                                          to=max(mat_graph_simi), 
                                          palette = 'Purples 3', rev = T),
             width = unit(1.5*ncol(mat_graph_simi)/nrow(mat_graph_simi), 'inch'),  
             height = unit(1.5, 'inch'), 
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (mat_graph_simi[i, j] > 0) {
                 grid.text(sprintf("%.02f", mat_graph_simi[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
               } 
             },
             use_raster=T, raster_by_magick = TRUE)
draw(p)
dev.off()



## Across archetypes
mat_graph_simi <- sapply(graph_combo_archetype, function(a) {
  sapply(graph_combo_archetype, function(b){
    jaccard_edgeset_similarity(a, b)
  })
})
diag(mat_graph_simi) <- 0
mat_graph_simi[upper.tri(mat_graph_simi)] <- 0
pdf(file.path(dir_res, sprintf('heatmap_graph_similarity_within_archetypes_%s.pdf', value_type)), 
    width = 1.5*ncol(mat_graph_simi)/nrow(mat_graph_simi) + 1.5, 
    height = 1.5+1.5, 
    useDingbats = F)
p <- Heatmap(mat_graph_simi, name = 'jaccard', 
             row_names_side = 'left', column_names_side = 'top',
             cluster_rows = F, cluster_columns = F, 
             col = heatmap_color_fun_cont(from=0, 
                                          to=max(mat_graph_simi), 
                                          palette = 'Purples 3', rev = T),
             width = unit(1.5*ncol(mat_graph_simi)/nrow(mat_graph_simi), 'inch'),  
             height = unit(1.5, 'inch'), 
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (mat_graph_simi[i, j] > 0) {
                 grid.text(sprintf("%.02f", mat_graph_simi[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
               } 
             },
             use_raster=F, raster_by_magick = TRUE)
draw(p)
dev.off()
