##
## Helper functions related to ecotype
##
## Yun Yan
##


suppressPackageStartupMessages({
  library(fs); library(readr); library(tidyverse); library(ggpubr)
  theme_set(theme_pubr(legend = 'right'))
  library(ComplexHeatmap)
  library(igraph)
  library(scales)
  library(colorspace)
  library(forcats)
  library(BiocParallel)
  library(Matrix)
  library(matrixStats)
})

#------ colors ------
  

pal_celltypes <- c(
  'tumor' = '#B4609D',
  'Endo' = '#7A333F',
  'Peri' = '#D3634D',
  'Fibro' = '#D3B839', 
  'Mye' = '#C0D576',
  'T' = '#2D2D70',
  'B' = '#5E91BC'
)
pal_celltypes <- c(
  'tumor' = '#D52126', # red
  'Endo' = 'burlywood4',
  'Peri' = 'purple',
  'Fibro' = '#FEE52C', # yellow
  'Mye' = '#117733', # green
  'T' = '#2F8AC4', # blue
  'B' = 'chocolate' # light blue
)
pal_pcr <- c(
  'pCR' = 'seagreen3',
  'RD'='tomato1', 
  'Unknown'='black', 'Excluded'='lightgrey')
pal_nmf4 <- c(
  'fNMF1' = '#0066CC',
  'fNMF2' = '#99cc00',
  'fNMF3' = '#ff9933',
  'fNMF4' = '#ff00cc'
)
pal_archetypes <- structure(as.character(pal_nmf4), names=paste0('ARC', 1:4))



cor_test_pvals_exclude_diag <- function(x, p.adjust.method='fdr', ...) {
  ## Like cor(), but gives the pvalues
  ## ...: to feed cor.test()
  ## The adjusted pvals will achieve the same result as psych::corr.test and lsr::correlate
  res_names <- colnames(x)
  n <- ncol(x)
  tri <- t(combn(1:n, m=2)) # upper half of matrix
  # tri <- rbind(tri, cbind(1:n, 1:n))    # diag indices
  val <- BiocParallel::bplapply(1:nrow(tri), function(i) {
    obj <- cor.test(x[, tri[i, 1]], x[, tri[i, 2]], ...)
    return(obj$p.value)
  })
  val <- c(unlist(val))
  if (!is.null(p.adjust.method)) {
    val <- p.adjust(p = val, method = p.adjust.method)
  }
  res <- sparseMatrix(i=tri[, 1], j=tri[, 2], x=val, dims=c(n, n))
  
  res <- as.matrix(res)
  
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  diag(res) <- 0
  rownames(res) <- colnames(res) <- res_names
  return(res)
}

co_geometric_mean <- function(x, na.rm=T) {
  ## Like cor(), but calculate the geometric mean
  res_names <- colnames(x)
  col_mean <- matrixStats::colMeans2(x, na.rm = na.rm)
  n <- ncol(x)
  tri <- t(combn(1:n, m=2))
  val <- apply(tri, 1, function(coord){
    # geometric mean of two values -- always cannot allow NAs
    # c(1, NA) = NA rather than 1
    exp(mean(log(c(col_mean[coord[1]], col_mean[coord[2]])), na.rm=FALSE))
  })
  res <- sparseMatrix(i=tri[, 1], j=tri[, 2], x=val, dims=c(n, n))
  res <- as.matrix(res)
  diag(res) <- col_mean
  
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  rownames(res) <- colnames(res) <- res_names
  return(res)
}
co_geometric_mean.one_obs <- function(x, na.rm=T) {
  # Like co_geometric_mean but runs on a single vector
  res_names <- names(x)
  col_mean <- as.numeric(x) #matrixStats::colMeans2(x, na.rm = na.rm)
  n <- length(x)
  tri <- t(combn(1:n, m=2))
  val <- apply(tri, 1, function(coord){
    # geometric mean of two values -- always cannot allow NAs
    # c(1, NA) = NA rather than 1
    exp(mean(log(c(col_mean[coord[1]], col_mean[coord[2]])), na.rm=FALSE))
  })
  res <- sparseMatrix(i=tri[, 1], j=tri[, 2], x=val, dims=c(n, n))
  res <- as.matrix(res)
  diag(res) <- col_mean
  
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  rownames(res) <- colnames(res) <- res_names
  return(res)
}

draw_thickness_legend <- function(actual_min, actual_max, viz_min, viz_max, ...){
  a <- seq(from=actual_min, to=actual_max, length.out=5)
  v <- rescale(a, to=c(viz_min, viz_max), from=c(actual_min, actual_max))
  legend(legend = sprintf('%.3f', a),
         lwd = v, ...
  )   
}

consensus_cluster_louvain_diagnosis <- function(
    graph, weights = NULL, resolution = 1, runs=10, verbose = F, ...) {
  # 
  # Values:
  # - a co-incidence matrix
  # - a vector of number of clusters
  
  node_names <- V(graph)$name
  mat <- matrix(0, nrow=length(node_names), ncol=length(node_names))
  rownames(mat) <- colnames(mat) <- node_names
  
  # diag(mat) <- runs
  num_clusters <- rep(0, length=runs)
  for (r in seq_len(runs)) {
    igraph_cluster_obj <- igraph::cluster_louvain(
      graph = graph, weights = weights, resolution = resolution)
    memb <- igraph::membership(igraph_cluster_obj)
    num_clusters[r] <- length(unique(memb))
    if (verbose) {print(table(memb))}
    for (memb_cat in unique(memb)) {
      nodes_in_cat <- names(memb)[memb == memb_cat]
      if (length(nodes_in_cat) == 1) {
        nodes_pair <- matrix(c(nodes_in_cat, nodes_in_cat), nrow = 1, ncol=2, byrow = T)
      } else {
        nodes_pair <- t(combn(nodes_in_cat, m=2))
      }
      for (rr in seq_len(nrow(nodes_pair))) {
        mat_i <- match(nodes_pair[rr, 1], node_names)
        mat_j <- match(nodes_pair[rr, 2], node_names)
        mat[mat_i, mat_j] <- mat[mat_i, mat_j] + 1
      }
    }
    
  }
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  diag(mat) <- runs ## force the feature itself is always self-coincident to help visualization
  return(list(mat = mat, 
              k = num_clusters))
}


consensus_cluster_louvain_member <- function(
    mat, igraph_cluster_method, ...){
  ig <- graph_from_adjacency_matrix(
    mat, mode = 'undirected', 
    weighted = TRUE, 
    diag = F
  )
  
  res <- igraph::membership(
    igraph_cluster_method(ig, ...))
  return(res)
}

weight.community=function(row,membership,weigth.within,weigth.between){
  # ref: https://stackoverflow.com/a/29098951
  if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
    weight=weigth.within
  }else{
    weight=weigth.between
  }
  return(weight)
}

zscore_na_rm <- function(x) {
  (x - mean(x, na.rm=T)) / sd(x, na.rm = T)
}

co_geometric_mean <- function(x, na.rm=T) {
  ## Like cor(), but calculate the geometric mean
  res_names <- colnames(x)
  col_mean <- matrixStats::colMeans2(x, na.rm = na.rm)
  n <- ncol(x)
  tri <- t(combn(1:n, m=2))
  val <- apply(tri, 1, function(coord){
    # geometric mean of two values -- always cannot allow NAs
    # c(1, NA) = NA rather than 1
    exp(mean(log(c(col_mean[coord[1]], col_mean[coord[2]])), na.rm=FALSE))
  })
  res <- sparseMatrix(i=tri[, 1], j=tri[, 2], x=val, dims=c(n, n))
  res <- as.matrix(res)
  diag(res) <- col_mean
  
  res[lower.tri(res)] <- t(res)[lower.tri(res)]
  rownames(res) <- colnames(res) <- res_names
  return(res)
}