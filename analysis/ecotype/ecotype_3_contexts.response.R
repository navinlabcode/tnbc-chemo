## Compare ecotype across patient groups (here: response)
##
## Yun Yan
##



dict_pat_to_response <- read.csv('./inputs/patient_metadata.csv')
dict_pat_to_response <- dict_pat_to_response[, c('patient', 'PCR_status')]


sample_names <- rownames(mat)

dict_sample <- dict_pat_to_response[sample_names, ]
dict_sample$sample_name <- dict_sample$patient
dict_sample$response <- dict_sample$PCR_status

#------------------ ~~~ start ~~~ --------------------
mat_use <- mat
rownames(mat_use)
mat_use_scaled <- apply(mat_use, 2, zscore_na_rm)
mat_use_value <- mat_use_scaled
cop_fea_run <- co_geometric_mean(mat_use_value)
cop_fea_run[is.na(cop_fea_run)] <- 0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Investigate response (PCR/nonPCR)            ~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------ separate samples to each context ------

dict_combo_sample_names <- dict_sample %>% dplyr::mutate(combo_id = paste0(response)) %>%
  dplyr::select(sample_name, combo_id) %>% 
  split(., f=.$combo_id)
dict_combo_sample_names <- lapply(dict_combo_sample_names, function(x) x$sample_name)
print(sapply(dict_combo_sample_names, length))
combo_names <- names(dict_combo_sample_names)
combo <- vector(mode = 'list', length = length(combo_names))
names(combo) <- combo_names

#------ creating cop for each context ------
for (zz in combo_names) {
  message(zz)
  sample_run <- dict_combo_sample_names[[zz]]
  
  print(length(sample_run))
  cop_fea_run <- co_geometric_mean(mat_use_value[sample_run, , drop=F])
  cop_fea_run[is.na(cop_fea_run)] <- 0
  cop_mat_use <- cop_fea_run
  # print(quantile(cop_mat_use, na.rm=T))
  Gf <- graph_from_adjacency_matrix(
    adjmatrix = cop_mat_use, mode = 'undirected', 
    weighted = TRUE, 
    diag = FALSE)
  edge_attr(Gf, 'label') <- apply(as_edgelist(Gf), 1, function(v) paste(v, collapse = '--'))
  
  print(quantile(E(Gf)$weight, na.rm = T))
  combo[[zz]] <- list(n_sample = length(sample_run),
                      cop_fea_matrix = cop_fea_run, 
                      graph = Gf)
}

combo_response <- combo

if (T) {
  legend_feature_ig_viz_weight_min <- min(
    sapply(combo, function(o) min(E(o$graph)$weight[E(o$graph)$label %in% igraph_feature_vizsig_edges])))
  legend_feature_ig_viz_weight_max <- max(
    sapply(combo, function(o) max(E(o$graph)$weight[E(o$graph)$label %in% igraph_feature_vizsig_edges])))
}


cat(legend_feature_ig_viz_weight_min, legend_feature_ig_viz_weight_max)

pdf(file.path(
  './rfiles/', 
  'igraph_contexts.%03d.pdf'),
  width = 7, height = 7, useDingbats = F, onefile = F)
for (zz in combo_names) {
  message(zz)
  cop_mat_use <- combo[[zz]]$cop_fea_matrix
  n_sample <- combo[[zz]]$n_sample
  Gf <- combo[[zz]]$graph
  print(quantile(E(Gf)$weight, na.rm=T))
  
  feature_ig_viz_weight_min = legend_feature_ig_viz_weight_min
  feature_ig_viz_weight_max = legend_feature_ig_viz_weight_max
  
  cat(feature_ig_viz_weight_min, feature_ig_viz_weight_max, '\n')
  
  # bad_edge <- (E(Gf)$weight<feature_ig_viz_weight_min)
  bad_edge <- (! E(Gf)$label %in% igraph_feature_vizsig_edges)
  bad_edge[is.na(bad_edge)] <- T
  
  Gf_viz <- delete_edges(Gf, E(Gf)[bad_edge])
  
  viz_edge_wt <- E(Gf_viz)$weight
  viz_edge_wt <- rescale(
    viz_edge_wt, to=c(feature_ig_plot_width_min, feature_ig_plot_width_max),
    from=c(legend_feature_ig_viz_weight_min, legend_feature_ig_viz_weight_max) )
  
  
  V(Gf_viz)$color <- pal_memb_feature[memb_feature[V(Gf_viz)]]#as.character(pal_celltypes[deframe(df_feature)[V(Gf_viz)$name]])
  # E(Gf_viz)$color <- 'grey'
  plot(
    Gf_viz,
    edge.color='slategray2', 
    edge.label=NA,
    vertex.label.color='black', #pal_celltypes[deframe(df_feature)[V(Gf_viz)$name]],
    vertex.size=4, vertex.label.cex=.6, vertex.label.dist = 0.4,
    edge.width = viz_edge_wt,
    edge.curved =F,
    layout=igraph_layout_feature,
    main = sprintf('%s (N=%d, top %s/%s edges)', 
                   zz, n_sample, length(E(Gf_viz)), length(E(Gf))) )
  legend("bottomright", legend=names(pal_memb_feature), col = pal_memb_feature,
         bty = "n", pch=20, cex=0.5, ncol=1,
         text.col='black', horiz = FALSE)
  draw_thickness_legend(
    legend_feature_ig_viz_weight_min, legend_feature_ig_viz_weight_max,
    feature_ig_plot_width_min, feature_ig_plot_width_max, 
    x = 'topleft', cex = .5, bty='n', title = 'co-presence')
  
}
dev.off()




save.image("./rfiles/ecotype.RData")
