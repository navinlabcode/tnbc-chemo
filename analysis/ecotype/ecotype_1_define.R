## Define Ecotypes
##
## Yun Yan
##


feature_ig_weight_cutoff <- 0

#------------------ ~~~ start ~~~ --------------------
corr_mat_use <- corr_fea
#------ igraph with the significance ------
Gf_qval <- graph_from_adjacency_matrix(
  adjmatrix = as.matrix(as.dist(corr_fea_pval)), mode = 'undirected', 
  weighted = TRUE, 
  diag = F)
edge_attr(Gf_qval, 'label') <- apply(as_edgelist(Gf_qval), 1, function(v) paste(v, collapse = '--'))
significant_edge_labels <- E(Gf_qval)$label[E(Gf_qval)$weight < 0.05]; str(significant_edge_labels)


#------ igraph with the correlation ------
Gf <- graph_from_adjacency_matrix(
  adjmatrix = as.matrix(as.dist(corr_mat_use)), mode = 'undirected', 
  weighted = TRUE, 
  diag = F)
edge_attr(Gf, 'label') <- apply(as_edgelist(Gf), 1, function(v) paste(v, collapse = '--'))

#------ determine the proper clustering resolution ------
Gf_cluster <- delete_edges(Gf, E(Gf)[weight<=feature_ig_weight_cutoff])
feature_ig_louvain_res_opts <- c(seq(from=0.01, to=10, length.out=1000))
feature_ig_louvain_res_opts <- sort(unique(feature_ig_louvain_res_opts))
set.seed(42)
o <- sapply(feature_ig_louvain_res_opts, function(feature_ig_louvain_res){
  igraph_cluster_obj <- igraph::cluster_louvain(
    Gf_cluster, resolution=feature_ig_louvain_res)
  memb_once <- igraph::membership(igraph_cluster_obj)
  tab_memb_once <- table(memb_once)
  
  return(c(feature_ig_louvain_res, length(unique(memb_once)), 
           max(tab_memb_once), min(tab_memb_once), median(tab_memb_once)))
} )
o <- t(o); head(o); o <- as.data.frame(o)
colnames(o) <- c('resolution', 'num_ecohub', 'max_ecohub_size', 'min_ecohub_size', 'med_ecohub_size')

y_metric <- 'num_ecohub'
y_metric_opts <- c('num_ecohub', 'max_ecohub_size', 'min_ecohub_size', 'med_ecohub_size')
o2 <- o[, c('resolution', 'num_ecohub', 'med_ecohub_size')] %>%
  pivot_longer(cols = c('num_ecohub', 'med_ecohub_size'), 
               names_to = 'metric', 
               values_to = 'n')
p <- ggplot(o2, aes(x=resolution, y=n)) + 
  geom_line(aes(color = metric))
print(p)

#------ choosing the 'best' resolution ------

feature_ig_louvain_res = 1.8 ## !!! MUST CHANGE ##
print(p + geom_vline(xintercept = feature_ig_louvain_res))

#------ defining ecotypes ------
set.seed(42)
igraph_cluster_diagnose <- consensus_cluster_louvain_diagnosis(
  graph = Gf_cluster, resolution = feature_ig_louvain_res, 
  run=1000
)
memb <- consensus_cluster_louvain_member(
  igraph_cluster_diagnose$mat, 
  igraph_cluster_method = igraph::infomap.community)
pal_memb <- try(structure(igraph::categorical_pal(n=length(unique(memb))), names=unique(memb)))

write_csv(enframe(memb), './rfiles/membership_ecotypes.csv')

if (T) {
  heatmap_clustering_method <- 'average'
  tmp <- factor(c(memb))
  pdf(file.path('./rfiles',
                'igraph.feature_correlations.diagnosis_consensus_clustering.pdf'),
      width = 8, height = 8, useDingbats = F, onefile = F)
  draw(
    Heatmap(igraph_cluster_diagnose$mat,
            # column_split  = 8, row_split = 8,
            clustering_method_columns = heatmap_clustering_method,
            clustering_method_rows = heatmap_clustering_method,
            right_annotation = rowAnnotation(df=data.frame(memb=tmp), 
                                             col=list(memb=pal_memb)), 
            width = unit(5, 'inch'), 
            height = unit(5, 'inch'), use_raster = F)
  )
  dev.off()
}

save.image("./rfiles/ecotype.RData")
