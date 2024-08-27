## Determine the layout of visualizing ecotypes
##
## Yun Yan
##


memb_feature <- memb
pal_memb_feature <- pal_memb


#------ adjust the edge width to fit the figure size ------
feature_ig_plot_width_min <- 0.1 
feature_ig_plot_width_max <- 3 


#------ Viz ecotypes using igraph ------
Gf_vizsig <- delete_edges(Gf_cluster, E(Gf_cluster)[!label %in% significant_edge_labels])
feature_ig_viz_weight_cutoff_min = quantile(E(Gf_vizsig)$weight, 0)
feature_ig_viz_weight_cutoff_max = quantile(E(Gf_vizsig)$weight, 1)
viz_edge_wt <- scales::rescale(
  E(Gf_vizsig)$weight, to=c(feature_ig_plot_width_min, feature_ig_plot_width_max), 
  from=c(feature_ig_viz_weight_cutoff_min, feature_ig_viz_weight_cutoff_max) )

V(Gf_vizsig)$color <- as.character(pal_memb_feature[memb_feature])#as.character(pal_celltypes[deframe(df_feature)[V(Gf_viz)$name]])
E(Gf_vizsig)$color <- 'maroon'
igraph_feature_vizsig_edges <- E(Gf_vizsig)$label

#------ Find the 'best' layout ------
feature_ig_viz_community_max_power = 18
pdf(file.path(
  './rfiles/', 
  'find_best_layout.%03d.pdf'),
  width = 7, height = 7, useDingbats = F, onefile = F)
for (feature_ig_viz_community_max_power in c(5, 10, 20, 50, 100, 18)) {
  set.seed(42)
  l_wt <- apply(
    get.edgelist(Gf_cluster), 1, 
    weight.community, memb_feature, feature_ig_viz_community_max_power, 1)
  l <- layout_with_fr(Gf_cluster,weights=l_wt)
  plot(
    Gf_vizsig, 
    # edge.color='grey', 
    edge.label=NA,
    vertex.label.color='black',#pal_memb[memb],
    vertex.size=4, 
    vertex.label.cex=0.6, 
    vertex.label.dist = 0.4,
    edge.width = viz_edge_wt,
    layout=l, 
    main = sprintf('graph clusters (%s edges)', length(E(Gf_vizsig)) ))
  legend("bottomright", legend=names(pal_memb_feature), col = pal_memb_feature,
         bty = "n", pch=20, cex=0.5, ncol=1,
         text.col='black', horiz = FALSE)
  draw_thickness_legend(
    feature_ig_viz_weight_cutoff_min, feature_ig_viz_weight_cutoff_max,
    feature_ig_plot_width_min, feature_ig_plot_width_max,
    x = 'topleft', cex = 1, bty='n', title = 'correlation')
}
dev.off()

#------ Using the 'best' layout ------
igraph_layout_feature <- l
V(Gf_vizsig)$color <- as.character(pal_memb_feature[memb_feature])#as.character(pal_celltypes[deframe(df_feature)[V(Gf_viz)$name]])
E(Gf_vizsig)$color <- 'lightsteelblue3'
# pdf(file.path('./rfiles/',
#               'igraph.feature_wise.pdf'), 
#     width = 8, height = 8, useDingbats = F, onefile = F)
plot(
  Gf_vizsig, 
  # edge.color='grey', 
  edge.label=NA,
  vertex.label.color='black',#pal_memb[memb],
  vertex.size=4, 
  vertex.label.cex=0.6, 
  vertex.label.dist = 0.4,
  edge.width = viz_edge_wt,
  layout=igraph_layout_feature, 
  main = sprintf('graph clusters (%s edges)', length(E(Gf_vizsig)) ))
legend("bottomright", legend=names(pal_memb_feature), col = pal_memb_feature,
       bty = "n", pch=20, cex=0.5, ncol=1,
       text.col='black', horiz = FALSE)
draw_thickness_legend(
  feature_ig_viz_weight_cutoff_min, feature_ig_viz_weight_cutoff_max,
  feature_ig_plot_width_min, feature_ig_plot_width_max,
  x = 'topleft', cex = 1, bty='n', title = 'correlation')
# dev.off()

save.image("./rfiles/ecotype.RData")
