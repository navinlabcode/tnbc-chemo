#
# Run AUCell on ST spots
# - archetype marker genes
# - metaprogram genes
#
# Yiyun Lin
#


#run in R4.0.3 singularity

aucell_wrapper <- function (seurat_inp, gene_list, assay='integrated', slot='data', 
                            thresh_select=NULL, trunc=TRUE, ncores=15, ...) {
  data_temp <- GetAssayData(seurat_inp, assay=assay, slot=slot)
  cat('Building ranks...\n')
  data_rank <- AUCell_buildRankings(data_temp, nCores=ncores, plotStats=FALSE, verbose=TRUE)
  cat('Calculating AUC...\n')
  data_AUC <- AUCell_calcAUC(gene_list, data_rank,
                             aucMaxRank=ceiling(0.2*nrow(data_rank)), nCores = ncores)
  data_AUC_df_raw <- t(data_AUC@assays@data$AUC) %>% data.frame(check.names = F)
  outlist <- list()
  outlist$AUC_raw <- data_AUC_df_raw
  if (trunc) {
    cat('Finding thresholds...\n')
    data_AUC_df_trunc <- data_AUC_df_raw
    data_ass <- AUCell_exploreThresholds(data_AUC, nCores=ncores, plotHist=FALSE, ...)
    cat('Generating thresholds matarix...\n')
    if (is.null(thresh_select)) {
      data_ass_df <- rep(sapply(data_ass, function(x) x$aucThr$selected), each=dim(data_AUC_df_raw)[1]) %>% 
        matrix(., nrow=dim(data_AUC_df_raw)[1], ncol=dim(data_AUC_df_raw)[2]) %>% data.frame(check.names = F)
    } else {
      data_ass_df <- rep(sapply(data_ass, function(x) x$aucThr$thresholds[thresh_select, 'threshold']), each=dim(data_AUC_df_raw)[1]) %>%
        matrix(., nrow=dim(data_AUC_df_raw)[1], ncol=dim(data_AUC_df_raw)[2]) %>% data.frame(check.names = F)
    }
    data_AUC_df_trunc[data_AUC_df_raw <= data_ass_df] <- 0
    outlist$AUC_trunc <- data_AUC_df_trunc
  }
  return(outlist)
}


if(T){
  options(stringsAsFactors = F)
  source("/volumes/USR1/yiyun/Script/scRNA/scRNA_packages.R")
  library("SChartPack")
  library("akima")
  library(patchwork)
  library("randomForestSRC")
  library("packcircles")
  library("dplyr")
  library(readr)
  library(AUCell)
  library("magrittr")
  library("dbscan")
  library("pheatmap")
  library("spatstat")
  # library("SeuratData")
  library(Seurat)
  library("reshape2")
  library("visNetwork")
  library("shiny")
  library("plotly")
  library("viridis")
  library("RColorBrewer")
  # library("ConsensusClusterPlus")
  library("philentropy")
  library('SChartPack')
  library('readr') 
  library(ggplot2)
  library(patchwork)
  # library('hdf5r',lib.loc = '/mnt/USR3/minhu/anaconda3/envs/seurat4/lib/R/library/')
}

# object bundles under seurat 3.2.0 in singularity rstudio_v4.0.3
objlist = read_rds('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_seuratRDS/object_list.rds') 

#calculate AUC run in 8787 R4.1.2
markers_NMF = read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/deliver/tumor_cells_only_20220707/deliver.archetypes_markers.top.rds')
markers_metatraits = read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/deliver/tumor_cells_only_20220707/deliver.mm_markers.rds')
markers = c(markers_NMF, markers_metatraits)
write_rds(markers,'/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/NMF_metatraits_marker.rds')

AUCell_list = lapply(objlist, function(x){AUCell=aucell_wrapper(x, gene_list = markers, assay='SCT', slot='data', thresh_select = 'Global_k1',ncores=100)})
write_rds(AUCell_list,'/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/NMF_metatraits_AUCres_022124.rds')


# NMF & Metatraits ---------------------------------------- 
#read AUC result add to metadata
AUCell_list = read_rds('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/NMF_metatraits_AUCres_022124.rds')
objlist_up = lapply(names(AUCell_list),function(pt){
  AUCell= AUCell_list[[pt]]
  AUCell_df <- AUCell$AUC_raw
  dim(AUCell_df)
  AUCell_df[1:5,1:5] %>% print()
  tumor = objlist[[pt]]
  for (i in 1:ncol(AUCell_df)) {
    cell_i <- colnames(AUCell_df)[i]
    tumor <- AddMetaData(tumor, metadata=AUCell_df[, i], col.name=paste0(cell_i, '_AUCell'))
    tumor@meta.data[tumor@meta.data$confidence != 'High_purity',paste0(cell_i, '_AUCell')]=0
  }
  return(tumor)
})
names(objlist_up) = names(objlist)

# featureplot
if(T){
  #NMF
  if(T){
    plist1 = lapply(names(objlist_up)[1:4],
                    function(i){
                      hl = AUCell_list[[i]]$AUC_raw %>% select(paste0('fNMF',1:4)) %>% 
                        summarise_if(is.numeric,max) %>% unlist() %>% median
                      plist_i = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 0.9,max.cutoff = hl,
                                                   features = c(paste0('fNMF',1:4,'_AUCell')),alpha = c(0.1, 1),
                                                   combine = F)
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      plist_i = lapply(plist_i,
                                       function(x){
                                         x+scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                                       })
                      p1 = plot_grid(plotlist = plist_i,nrow = 2)
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    plist2 = lapply(names(objlist_up)[c(5:10,12)],
                    function(i){
                      hl = AUCell_list[[i]]$AUC_raw %>% select(paste0('fNMF',1:4)) %>% 
                        summarise_if(is.numeric,max) %>% unlist() %>% median
                      plist_i = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 4,max.cutoff = hl,
                                                   features = c(paste0('fNMF',1:4,'_AUCell')),alpha = c(0.1, 1),
                                                   combine = F)
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      plist_i = lapply(plist_i,
                                       function(x){
                                         x+scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                                       })
                      p1 = plot_grid(plotlist = plist_i,nrow = 1)
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    plist3 = lapply(names(objlist_up)[c(11)],
                    function(i){
                      hl = AUCell_list[[i]]$AUC_raw %>% select(paste0('fNMF',1:4)) %>% 
                        summarise_if(is.numeric,max) %>% unlist() %>% median
                      plist_i = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 2.5,max.cutoff = hl,
                                                   features = c(paste0('fNMF',1:4,'_AUCell')),alpha = c(0.1, 1),
                                                   combine = F)
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      plist_i = lapply(plist_i,
                                       function(x){
                                         x+scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                                       })
                      p1 = plot_grid(plotlist = plist_i,nrow = 1)
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    pdf('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/fNMF_summary/fNMF_ST_022124.pdf',width = 10,height = 8)
    plist1 %>% print()
    plist2 %>% print()
    plist3 %>% print()
    dev.off()
  }
  
  #metatraits
  if(T){
    plist1 = lapply(names(objlist_up)[1:4],
                    function(i){
                      metatraits = grep(pattern = 'M\\d{1,2}',colnames(AUCell_list[[i]]$AUC_raw),value = T)
                      hl = AUCell_list[[i]]$AUC_raw %>% 
                        select(all_of(metatraits)) %>% 
                        summarise_if(is.numeric,max) %>% unlist() %>% median
                      plist_i = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 0.9,max.cutoff = hl,
                                                   features = c(paste0(metatraits,'_AUCell')),alpha = c(0.1, 1),
                                                   combine = F)
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      plist_i = lapply(plist_i,
                                       function(x){
                                         x+scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                                       })
                      p1 = plot_grid(plotlist = plist_i,ncol = 4)
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    plist2 = lapply(names(objlist_up)[c(5:10,12)],
                    function(i){
                      metatraits = grep(pattern = 'M\\d{1,2}',colnames(AUCell_list[[i]]$AUC_raw),value = T)
                      hl = AUCell_list[[i]]$AUC_raw %>% 
                        select(all_of(metatraits)) %>% 
                        summarise_if(is.numeric,max) %>% unlist() %>% median
                      plist_i = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 4,max.cutoff = hl,
                                                   features = c(paste0(metatraits,'_AUCell')),alpha = c(0.1, 1),
                                                   combine = F)
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      plist_i = lapply(plist_i,
                                       function(x){
                                         x+scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                                       })
                      p1 = plot_grid(plotlist = plist_i,nrow = 2)
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    plist3 = lapply(names(objlist_up)[c(11)],
                    function(i){
                      metatraits = grep(pattern = 'M\\d{1,2}',colnames(AUCell_list[[i]]$AUC_raw),value = T)
                      hl = AUCell_list[[i]]$AUC_raw %>% 
                        select(all_of(metatraits)) %>% 
                        summarise_if(is.numeric,max) %>% unlist() %>% median
                      plist_i = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 2.5,max.cutoff = hl,
                                                   features = c(paste0(metatraits,'_AUCell')),alpha = c(0.1, 1),
                                                   combine = F)
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      plist_i = lapply(plist_i,
                                       function(x){
                                         x+scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                                       })
                      p1 = plot_grid(plotlist = plist_i,nrow = 2)
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    pdf('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/fNMF_summary/metatraits_ST.pdf',width = 10,height = 8)
    plist1%>% print()
    plist2%>% print()
    plist3%>% print()
    dev.off()
  }
}

# cancer genes ----------------------------------------
# add genes norm expr to metadata

genes = c('CCND1', 'MYC', 'MTDH', 'CDKN2A', 'CKS1B', 'GADD45A', 'GATA3')


objlist_up = lapply(names(objlist),function(pt){
  tumor = objlist[[pt]]
  mat = GetAssayData(tumor,assay = 'SCT',slot = 'data')
  df_mat = mat %>% t() %>%  data.frame() %>% select(all_of(genes))
  for (i in 1:ncol(df_mat)) {
    cell_i <- colnames(df_mat)[i]
    tumor <- AddMetaData(tumor, metadata=df_mat[, i], col.name=paste0(cell_i,'_SCT'))
    tumor@meta.data[tumor@meta.data$confidence != 'High_purity',paste0(cell_i,'_SCT')]=0
  }
  return(tumor)
})

names(objlist_up) = names(objlist)
gmat_list = lapply(names(objlist),function(pt){
  tumor = objlist_up[[pt]]
  genemat = tumor@meta.data %>% select(all_of(paste0(genes,'_SCT')))
  genemat$patient = pt
  return(genemat)
})
gmat = do.call(rbind,gmat_list)
hlall =gmat %>%  data.frame() %>% dplyr::group_by(patient) %>% 
  summarise_if(is.numeric,max) %>% select(-patient) %>%  
  data.frame() %>% apply(.,2,median)

#featureplot
if(T){
  for(g in names(hlall)){
    hl =hlall[g]
    plist1 = lapply(names(objlist_up)[1:4],
                    function(i){ 
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      p1 = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 0.9,max.cutoff = hl,
                                                   features = g,alpha = c(0.1, 1),slot = 'data',
                                                   combine = F)[[1]]+
                        scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                      
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    plist2 = lapply(names(objlist_up)[c(5:10,12)],
                    function(i){
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      p1 = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 4,max.cutoff = hl,
                                              features = g,alpha = c(0.1, 1),slot = 'data',
                                              combine = F)[[1]]+
                        scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
            
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    plist3 = lapply(names(objlist_up)[c(11)],
                    function(i){
                      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
                      p1 = SpatialFeaturePlot(objlist_up[[i]],pt.size.factor = 2.5,max.cutoff = hl,
                                              features = g,alpha = c(0.1, 1),slot = 'data',
                                              combine = F)[[1]]+
                        scale_fill_gradientn(colours = myPalette(100),limits = c(0,hl))
                     
                      title_plot <- ggdraw() + draw_label(paste0(unique(objlist_up[[i]]$patient),
                                                                 '_',unique(objlist_up[[i]]$pCR)))
                      p = plot_grid(title_plot, p1, ncol=1, rel_heights=c(0.1, 1))
                      return(p)
                    })
    
    pdf(paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/fNMF_summary/',g,'_expression.pdf'),width = 8,height = 8)
    plist1 %>% print()
    plist2 %>% print()
    plist3 %>% print()
    dev.off()
  }
  
}


