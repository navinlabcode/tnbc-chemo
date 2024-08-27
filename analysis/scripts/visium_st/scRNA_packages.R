# Helper functions
#
# Yiyun Lin
#

library(dplyr)
library(stringr)
library(janitor)
library(Seurat)
library(future)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape)
library(plyr)
library(CellChat)

CreIntegObject <- function(ipdir, sample, timepoint,tissue.type){
  print("Read10X expression matrix")
  ART44_Pre.data <- Read10X(data.dir = paste0(ipdir,"outs/filtered_feature_bc_matrix/"))
  print("Create Seurat object")
  ART44_Pre_C <- CreateSeuratObject(counts = ART44_Pre.data, project = paste0(sample), min.cells = 3, min.features = 200)
  ART44_Pre_C[["percent.mt"]] <- PercentageFeatureSet(ART44_Pre_C, pattern = "^MT-")
  ART44_Pre_C[["percent.rb"]] <- PercentageFeatureSet(ART44_Pre_C, pattern = "^RP")
  ART44_Pre_C[["percent.RBC"]] <- PercentageFeatureSet(ART44_Pre_C, pattern = "^HBB|^HBA1|^HBA2")
  p1<-VlnPlot(ART44_Pre_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.RBC"), ncol = 5)
  print(p1)
  pdf(file = paste0(sample,"QC.pdf"), width = 15, height = 7)
  print(p1)
  dev.off()
 
  sample.name = stringr::str_extract(string = sample, pattern = '^ART\\d{2,3}');sample.name
  tp = timepoint
  tissue.type = tissue.type
  metaname = paste0(sample.name,tp,tissue.type)
  
  print(paste0("Set object name as ", metaname))
  ART44_Pre_C@meta.data$sample <- sample
  ART44_Pre_C@meta.data$timepoint <- timepoint
  print("Normalize Data...")
  ART44_Pre_C <- NormalizeData(ART44_Pre_C,normalization.method = "LogNormalize", scale.factor = 10000)
  print("Scale Data...")
  print("Find variable features by variance stabilizing transformation (vst), before finding anchors")
  ART44_Pre_C <- FindVariableFeatures(ART44_Pre_C, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(ART44_Pre_C)
  print("use all gene for scaling")
  ART44_Pre_C <- ScaleData(ART44_Pre_C, features = all.genes, verbose = F)
  print("Running PCA, PCs 1-30")
  ART44_Pre_C <- RunPCA(ART44_Pre_C, features = VariableFeatures(object = ART44_Pre_C))
  print("FindNeighbors for samples")
  ART44_Pre_C <- FindNeighbors(ART44_Pre_C, dims = 1:30)
  print("FindClusters for samples")
  ART44_Pre_C <- FindClusters(ART44_Pre_C, resolution = 0.5)
  print("Running TSNE and UMAP, PCs 1-30")
  # ART44_Pre_C <- RunTSNE(ART44_Pre_C, dims = 1:30)
  ART44_Pre_C <- RunUMAP(ART44_Pre_C, dims = 1:30)
  ART44_Pre_C <- RenameCells(ART44_Pre_C, new.names = paste0(timepoint,"_",colnames(x = ART44_Pre_C)))
  return(ART44_Pre_C)
}

#select PCs dims -----
select_dim = function(obj){
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1   
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  co2
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  pcs %>% return()
}

### CopyKat ###
# run command in terminal #
# Rscript /volumes/lab/users/ruligao/test_COPYCAT/bin/run_COPYCAT_main_V17.R --path.V3 /volumes/seq/projects/ARTEMIS/cellranger_v3_10/ART44-PreTx-Core/outs/filtered_feature_bc_matrix/ --cores 10
#allign new identity
PredTum <- function(predir, name, tp, ART, ifmanual){
  if(ifmanual==F){
    pred_results <- read.table(paste0(predir,name,"/copycat/KS.cut=0.05/",name,"_prediction.txt"))
    pred_results = pred_results[-1,]
    pred <- as.data.frame(pred_results[,2])
    if(sum(stringr::str_detect(pred_results[,1],tp))<nrow(pred_results)){
      rownames(pred)[which(stringr::str_detect(pred_results[,1],tp,negate = T))] <- paste0(tp,'_', pred_results[,1])
      rownames(pred)[which(stringr::str_detect(pred_results[,1],tp))] <- pred_results[,1]
      }else{rownames(pred) = pred_results[,1]}
    sID <- rep(name, nrow(pred))
    pred <- cbind(pred,sID)
    dim(pred); colnames(pred) <- c("class","sample.name")
    pred$class = gsub(pattern = 'c2:low.confidence',replacement =  'nondiploid', pred$class)
    pred$sample.name = gsub(pattern = 'c1:low.confidence',replacement =  'diploid', pred$sample.name)
    ART <- AddMetaData(ART, pred)
    ART@meta.data$class = factor(ART@meta.data$class, levels =c("nondiploid","diploid","not_predicted"))
    ART@meta.data$class[which(is.na(ART@meta.data$class))] <- factor("not_predicted")
    print(table(ART@meta.data$class))
    return(ART)
  }else{
    sample.name = stringr::str_extract(string = name, pattern = '^ART\\d{2,3}')
    pred_results <- read.table(paste0('/volumes/lab/users/yiyun/Project/ARTEMIS_RNA/',sample.name,'/',sample.name,"_predictionM.txt"))
    rownames(pred_results) = gsub('MidTx_MidTx_', 'MidTx_',rownames(pred_results))
    rownames(pred_results) = gsub('PreTx_PreTx_', 'PreTx_',rownames(pred_results))
    rownames(pred_results) = gsub('PostTx_PostTx_', 'PostTx_',rownames(pred_results))
    pred_results$cellname =rownames(pred_results)
    pred_results
    pred <- as.data.frame(pred_results[,c('group','cluster')])
    rownames(pred) <- rownames(pred_results)
    #rownames(pred) <- pred_results$V1
    colnames(pred) <- c("class", "cluster")
    sID <- rep(name, nrow(pred))
    pred <- cbind(pred,sID)
    pred$cluster <- paste0(sample.name,"_", pred$cluster)
    dim(pred)
    #pred <- pred[-1,] #remove ruli's header
    ART <- AddMetaData(ART, pred)
    levels(ART@meta.data$class) <- factor(c("Aneuploid","Diploid","Not_predicted"))
    ART@meta.data$class[is.na(ART@meta.data$class)] = factor("Not_predicted")
    ART@meta.data$cluster[is.na(ART@meta.data$cluster)] <- paste0(sample.name,"_NP")
    print(table(which(ART@meta.data$cluster==1)))
    print(table(ART@meta.data$class))
    print(head(pred))
    return(ART)
  }
}



str_clean <- function(x){
  name=x
  name = gsub('\\/','_',name)
  name = gsub('\\+','_',name)
  name = gsub('\\ ','',name)
  name = gsub('__','_',name)
  name = gsub('-','_',name)
  name = gsub('\\(','_',name)
  name = gsub('\\)','_',name)
  return(name)
}


timepoint <-function(x){
  if(length(x)==1){
    PreTx = stringr::str_detect(x, regex("pre", ignore_case = TRUE))
    MidTx = stringr::str_detect(x, regex("mid", ignore_case = TRUE))
    PostTx = stringr::str_detect(x, regex("post", ignore_case = TRUE))
    if(PreTx){timepoint = "PreTx";}else{if(MidTx){timepoint = "MidTx"}else{timepoint = "PostTx"}}
    return(timepoint)
  }else{
    tp = function(x){
      PreTx = stringr::str_detect(x, regex("pre", ignore_case = TRUE))
      MidTx = stringr::str_detect(x, regex("mid", ignore_case = TRUE))
      PostTx = stringr::str_detect(x, regex("post", ignore_case = TRUE))
      if(PreTx){timepoint = "PreTx";}else{if(MidTx){timepoint = "MidTx"}else{timepoint = "PostTx"}}
      return(timepoint)}
    timepoint = c()
    for(i in 1:length(x)){
      b = x[i]
      a = tp(b)
      timepoint[i] =a
    }
    return(timepoint)
  }
}

tissue <-function(x){
  Core = stringr::str_detect(x, regex("_C|-C|core", ignore_case = TRUE))
  FNA = stringr::str_detect(x, regex("_F|-F|FNA", ignore_case = TRUE))
  Sx = stringr::str_detect(x, regex("_Sx|TxSx|Sur|ART75_Post", ignore_case = TRUE))
  if(Core){tissue = "Core"}else{if(FNA){tissue = "FNA"}else{if(Sx){tissue = "Sx"}else{tissue = "Core"}}}
  return(tissue)
}


reorderTP <- function(x){
  a = unique((x))
  b = c("PreTx","MidTx","PostTx")
  c=c()
  for(i in b){
    c = c(c,grep(i,a,value = T))
  }
  message(paste0("sequence reordered as "))
  print(c)
  return(c)
}

########################################  aucell_wrapper  ##############################################  
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


library(RColorBrewer)
########################################  colormaker  ##############################################  
colormaker <- function(x, pattern = "Set2"){
  n_color = length(unique(x)); 
  if(n_color>8){
    minor_col<-setNames(colorRampPalette(brewer.pal(8, pattern))(n_color),sort(unique(x),decreasing = F)) 
  }else{minor_col<-setNames(brewer.pal(n_color, pattern),sort(unique(x),decreasing = F))}
  return(minor_col[1:n_color])
}

tpcolorset <- function(x){
  Groups_set <- c("PreTx", "MidTx", "PostTx")
  col<-c("#60c0c0", "orange", "#f07878")
  col = as.vector(col)
  names(col) <- Groups_set
  
  a = levels(x) %>% as.character()
  c= c()
  for(i in 1:length(a)){
    c=c(c,timepoint(a[i]))
  }
  c = col[names(col) %in% c]
  return(c)
}



Makersumary <- function(x,int.met,res = 0.3,step = 2){
  message(paste0('Default resolution is ',res))
  res.0.1 = as.data.frame(Idents(x)) %>% set_colnames(paste0("res.",res))
  x <- AddMetaData(x, res.0.1)
  GenMarker <- c("PTPRC", "CD3D", "CD8A", "CD4","GZMB","CD19", "CD14","CD68",'CD1C', 'FCER1A', "ACTA2",
                 "DCN","EPCAM","KRT18","KRT5","KRT14",'PECAM1',"PLVAP","LYVE1")
  dat =GetAssayData(x, slot = "data")
  mat=do.call(rbind,lapply(sort(unique(res.0.1[,1])), function(t){rowMeans(dat[,res.0.1[,1]==t])}))
  mat=scale(mat)
  fi = data.frame('Gene'=c('LUM',"DCN",'FAP',"COL1A1",'PDPN',"APOD","COL1A2","COL3A1",'COL6A1',"ACTA2",'MYLK','MYL9','ACTG2','MYH11'),'Cell.type' = 'Fibroblast')
  t = data.frame('Gene'=c('PTPRC','CD2','CD3D','CD3E','CD3G','CD8A','NKG7','CD4','FOXP3',"CTLA4"),'Cell.type' = 'Tcell')
  m_mp_d = data.frame('Gene'= c('CD14','LYZ','MAFB','CSF1R','MSR1','CD300E','CD1C'),'Cell.type' ='Mono/Mq/DC')
  pDC = data.frame('Gene'= c('LILRA4','IRF7','PACSIN1','TCF4','IL3RA','GZMB'),'Cell.type' = 'pDC')
  neu =  data.frame('Gene'= c('CXCL8','CXCR2','FCGR3B','CSF3R','S100A8','MMP8'),'Cell.type' ='Neutrophil')
  # m = data.frame('Gene'=c("CD14",'FCGR3A',"CD68","FCGR1A","S100A9","HLA-DRA",'VCAN','CSF1R'),'Cell.type' = 'Mono')
  # dc = data.frame('Gene'=c('CD1C', 'IL3RA','FCER1A','CLEC10A','FCGR2B'),'Cell.type' = 'DC')
  b = data.frame('Gene'=c('CD79A','CD79B','CD19','MS4A1','BANK1','FCRL5'),'Cell.type' = 'Bcell')
  naiveB =data.frame('Gene'=c('MS4A1','CXCR4','HLA-DRA','IL4R'), 'Cell.type' = 'naiveB')
  plsB =data.frame('Gene'=c("IGKC",'IGHA1',"IGHM",'JCHAIN','IGHG1'), 'Cell.type' = 'Plasma')
  endo = data.frame('Gene'=c("PECAM1",'GJA5','VWF',"PLVAP","LYVE1","HSPG2",'ENG'),'Cell.type' = 'Endo')
  mast = data.frame('Gene'=c("TPSB2","TPSAB1"),'Cell.type' = 'Mast')
  epi = data.frame('Gene'=c('EPCAM', 'EGFR', 'CDH1'),'Cell.type' = 'Epithelial')
  basal = data.frame('Gene'=c('KRT14','ITGA6','KRT5','TP63','KRT17'),'Cell.type' = 'Basal')
  lum = data.frame('Gene'=c('MME','KRT8','KRT18','KRT19','FOXA1','GATA3', 'MUC1','CD24'),'Cell.type' = 'Luminal')
  lumP = data.frame('Gene'=c('KIT','GABRP'),'Cell.type' = 'Lum_Pro')
  RBC =  data.frame('Gene'=c('HBB','HBA1', 'HBA2','HEMGN','ALAS2'),'Cell.type' = 'RBC')
  gene = rbind(t,m_mp_d,pDC,neu,b,plsB,endo,mast,fi,epi,basal,lum,lumP,RBC); rownames(gene)= gene$Gene
  gene= gene[gene$Gene%in%colnames(mat),]
  which(colnames(mat)%in%gene$Gene)
  gmat = mat[,which(colnames(mat)%in%gene$Gene)]
  gmat = gmat[,gene$Gene %>% as.character()]; rownames(gmat) = as.character(sort(unique(res.0.1[,1])))
  dim(gmat);dim(gene)
  library(ComplexHeatmap)
  anno_col = gene$Cell.type %>% data.frame() %>% set_rownames(gene$Gene)
  cols = colormaker(x$seurat_clusters)[1:length(unique(gene$Cell.type))];names(cols) = unique(gene$Cell.type)
  ha = HeatmapAnnotation(cellmarker = anno_col$.,col = list(cellmarker = cols),show_annotation_name = T)
  ht = Heatmap(gmat,cluster_rows = T,cluster_columns = F,
               clustering_distance_rows = "pearson",
               clustering_method_rows = 'complete',
               name = 'AvgExpression',column_names_rot = 45,
               column_names_side = "bottom",border = T,
               top_annotation =ha, 
               column_split = gene$Cell.type,
               col = circlize::colorRamp2(breaks =c(-2,0,2),c("dodgerblue3", "white", "firebrick3")))
  
  hm =grid.grabExpr(draw(ht))
  pe = plot.new()
  pdf(file = paste0(step,"_marker_gene_summary_","_",int.met,"_",res,".pdf"), width = 20, height = 20)
  #print(FeaturePlot(x, features = GenMarker,label =T,n=5))
  print(plot_grid(pe,hm, rel_widths = c(1,10),rel_height = 5))
  dev.off()
  print(plot_grid(pe,hm, rel_widths = c(1,10)))
  return(gmat)
}

subMarkersumary <- function(x,Tmarker,int.met = 'RNA',res = 0.3,step = 2){
  message(paste0('Default resolution is ',res))
  res.0.1 = as.data.frame(Idents(x)) %>% set_colnames(paste0("res.",res))
  DefaultAssay(x) = 'RNA'
  x <- AddMetaData(x, res.0.1)
  dat =GetAssayData(x, slot = "data")
  mat=do.call(rbind,lapply(sort(unique(res.0.1[,1])), function(t){rowMeans(dat[,res.0.1[,1]==t])}))
  #mat=scale(mat)
  gene = data.frame()
  for(i in 1:length(Tmarker)){
    df = data.frame(Gene = Tmarker[[i]],Cell.type = names(Tmarker)[i])
    gene = rbind(gene,df)
  }

  gene= gene[gene$Gene%in%colnames(mat),]
  which(colnames(mat)%in%gene$Gene)
  gmat = mat[,which(colnames(mat)%in%gene$Gene)]
  gmat = gmat[,gene$Gene %>% as.character()]; rownames(gmat) = as.character(sort(unique(res.0.1[,1])))
  dim(gmat);dim(gene)
  library(ComplexHeatmap)
  
  cols = colormaker(x$seurat_clusters)[1:length(unique(gene$Cell.type))];names(cols) = unique(gene$Cell.type)
  ha = HeatmapAnnotation(cellmarker = gene$Cell.type,col = list(cellmarker = cols),show_annotation_name = T)
  ht = Heatmap(scale(gmat),cluster_rows = T,cluster_columns = F,
               clustering_distance_rows = "pearson",
               clustering_method_rows = 'complete',
               name = 'AvgExpression',column_names_rot = 45,
               column_names_side = "bottom",border = T,
               top_annotation =ha, 
               column_split = gene$Cell.type,
               col = circlize::colorRamp2(breaks =c(-2,0,2),c("dodgerblue3", "white", "firebrick3")))
  
  hm =grid.grabExpr(draw(ht))
  pe = plot.new()
  pdf(file = paste0(step,"_marker_gene_summary_",res,".pdf"), width = 20, height = 20)
  #print(FeaturePlot(x, features = GenMarker,label =T,n=5))
  print(plot_grid(pe,hm, rel_widths = c(1,10),rel_height = 5))
  dev.off()
  print(plot_grid(pe,hm, rel_widths = c(1,10)))
  return(gmat)
}

subMarkersumary_high <- function(x,Tmarker,split=4,int.met = 'RNA',meta=NA){
  res.0.1 = as.data.frame(Idents(x)) %>% set_colnames(paste0("res."))
  DefaultAssay(x) = int.met
  x <- AddMetaData(x, res.0.1)
  dat =GetAssayData(x, slot = "data")
  mat=do.call(rbind,lapply(sort(unique(res.0.1[,1])), function(t){rowMeans(dat[,res.0.1[,1]==t])}))
  #mat=scale(mat)
  gene = data.frame()
  for(i in 1:length(Tmarker)){
    df = data.frame(Gene = Tmarker[[i]],Cell.type = names(Tmarker)[i])
    gene = rbind(gene,df)
  }
  
  gene= gene[gene$Gene%in%colnames(mat),]
  which(colnames(mat)%in%gene$Gene)
  gmat = mat[,which(colnames(mat)%in%gene$Gene)]
  gmat = gmat[,gene$Gene %>% as.character()]; rownames(gmat) = as.character(sort(unique(res.0.1[,1])))
  dim(gmat);dim(gene)
  library(ComplexHeatmap)

  cols = colormaker(gene$Cell.type)[1:length(unique(gene$Cell.type))];names(cols) = unique(gene$Cell.type)
  ha = rowAnnotation(a = gene$Cell.type,col = list(a = cols),show_annotation_name = T)
  data = t(gmat)
  if(split==0){
    ht = Heatmap(t(scale(t(data))),cluster_rows = F,cluster_columns = F,
                 # clustering_distance_columns  = "pearson",
                 # clustering_method_columns  = 'complete',
                 name = 'Scaled Expression',column_names_rot = 45,
                 column_names_side = "bottom",border = F,
                 left_annotation =ha, 
                 row_split = gene$Cell.type,
                 column_split = meta,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 col = circlize::colorRamp2(breaks =c(-2,0,2),c("dodgerblue3", "white", "firebrick3")))
  }else{
    ht = Heatmap(t(scale(t(data))),cluster_rows = F,cluster_columns = T,
                 clustering_distance_columns  = "pearson",
                 clustering_method_columns  = 'complete',
                 name = 'Scaled Expression',column_names_rot = 45,
                 column_names_side = "bottom",border = F,
                 left_annotation =ha, 
                 row_split = gene$Cell.type,
                 column_split = split,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 col = circlize::colorRamp2(breaks =c(-2,0,2),c("dodgerblue3", "white", "firebrick3")))
  }
  return(ht)
  ht
}


cell_freq_plot <- function(x,variable,levels,xasis,ifPlt = T,ifTME ='TME'){
  if(ifTME == 'TME'){
    df = x@active.ident
    levels2 = levels(x@active.ident)
  }
  if(ifTME == 'CNA'){
    df  = x$class
    levels2 = sort(unique(x$class))
  }
  if(ifTME == 'epi'){
    df = x$epiID
    levels2 = sort(unique(x$epiID))
  }
  x_freq<- data.frame(melt(table(df ,variable)))
  x_freq$Var.1 <- factor(x_freq$df , levels = levels2)
  x_freq$Var.2 <- factor(x_freq$variable, levels = levels)
  x_freq <- ddply(x_freq,.(Var.2),transform, percent=value/sum(value)*100)
  x_freq = ddply(x_freq, .(Var.2), transform, pos = (cumsum(value) - 0.5 * value))
  x_freq$label = paste0(sprintf("%.0f", x_freq$percent), "%")
  p1 <- ggplot(x_freq, aes(y = percent, x = Var.2, fill = Var.1)) + 
    geom_bar(position=position_stack(),stat="identity")+
    labs(x=xasis,y="Percentage", fill="Cell Types")+
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4)+
    theme(axis.text.x = element_text(size = 14, face = "bold", angle = 90,vjust = 0.5),axis.text.y = element_text(size = 15),axis.title=element_text(size=16,face = "bold"))+
    scale_fill_manual(name="Celltype",values = colormaker(unique(df)))

  if(ifPlt == T){return(p1)}
  if(ifPlt == F){return(x_freq)}
}

FindmarkerResponse = function(seu,timepoint,test ='DESeq2'){
  seu = subset(seu,subset = timepoint==timepoint)
  Idents(seu) = seu$response
  if(test == 'DESeq2'){
    seu[["RNA"]]@counts<-as.matrix(seu[["RNA"]]@counts)+1
    MidTxBmarker = FindMarkers(seu, ident.1 = "NonResponder", ident.2 = "Responder", test.use = "DESeq2")
  }else{
    MidTxBmarker = FindMarkers(seu, ident.1 = "NonResponder", ident.2 = "Responder")
  }
  MidTxBmarker = MidTxBmarker[order(MidTxBmarker$avg_logFC,decreasing = T),]
  return(MidTxBmarker)
}


addID= function(soj,ID){
  soj <- RenameIdents(soj, ID)
  table(Idents(soj))
  soj$ST.Ident1 = Idents(soj)
}


addRes = function(x){
  Response = rep(NA, length(x$patient)) %>%  magrittr::set_names(colnames(x))
  Response[which(x$patient %in% c('ART44','ART45','ART46','ART100'))] ='RCBIII'
  Response[which(x$patient %in% c('ART02','ART13','ART21','ART43','ART55','ART63','ART75','ART77','ART66','ART89','ART170','ART159','ART155'))] ='RCBII'
  Response[which(x$patient %in% c('ART10','ART17','ART51','ART70','ART92','ART94','ART134','ART180','ART136','ART153','ART142'))] = 'RCBI'
  Response[which(x$patient %in% c('ART05','ART07','ART12','ART19','ART23','ART30','ART31','ART50','ART53','ART60','ART62',
                                  'ART71','ART90','ART96','ART104','ART124','ART163','ART174','ART181','ART168','ART161','ART133','ART140','ART132','ART118','ART122','ART117'))] = 'PCR'
  Response[is.na(Response)] = 'Unknown'
  print(table(Response))
  
  GenResponse = rep(NA, length(x$patient)) %>% magrittr::set_names(colnames(x))
  GenResponse[which(x$patient %in% c('ART44','ART45','ART46','ART100',
                                     'ART02','ART13','ART21','ART43','ART55','ART63','ART75','ART77','ART66','ART89','ART170','ART159','ART155', 
                                     'ART10','ART17','ART51','ART70','ART92','ART94','ART134','ART180','ART136','ART153','ART142'))] ='NonResponder'
  GenResponse[which(x$patient %in% c('ART05','ART07','ART12','ART19','ART23','ART30','ART31','ART50','ART53','ART60','ART62','ART71','ART90',
                                     'ART96','ART104','ART124','ART163','ART174','ART181','ART168','ART161','ART133','ART140','ART132','ART118','ART122','ART117'))] = 'Responder'
  GenResponse[is.na(GenResponse)] = 'Unknown'
  print(table(GenResponse))
  
  Chem = rep(NA, length(x$patient)) %>% magrittr::set_names(colnames(x))
  Chem[which(x$patient %in% c('ART02','ART05','ART07','ART10','ART12','ART13','ART18','ART19','ART44','ART62'))] ='V2'
  Chem[is.na(Chem)] = 'V3'
  print(table(Chem))
  
  Vanderbilt = rep(NA, length(x$patient)) %>% magrittr::set_names(colnames(x))
  Vanderbilt[which(x$patient %in% c('ART17','ART18','ART44','ART43','ART50','ART70'))] ='BL1'
  Vanderbilt[which(x$patient %in% c('ART08','ART62','ART71','ART100','ART104',"ART157"))] ='BL2'
  Vanderbilt[which(x$patient %in% c('ART02','ART07','ART21','ART19','ART23','ART30','ART31',
                                    'ART39','ART53','ART60',"ART62",'ART78','ART90','ART96',
                                    'ART125','ART156','ART163','ART179','ART181',"ART192", 
                                    "ART142", "ART140", "ART114", "ART117", "ART190", "ART188","ART186"))] ='IM'
  Vanderbilt[which(x$patient %in% c('ART89','ART174',"ART119", "ART115","ART187"))] ='LAR'
  Vanderbilt[which(x$patient %in% c('ART05','ART12','ART13','ART46','ART51','ART55','ART66','ART75',
                                    'ART92','ART94','ART124','ART134','ART170',"ART194","ART159","ART153","ART132"))] ='M'
  Vanderbilt[which(x$patient %in% c("ART45","ART155"))] ='MSL'
  Vanderbilt[is.na(Vanderbilt)] = 'UNS'
  print(table(Vanderbilt))
  
  Res=cbind(Response,GenResponse,Chem,Vanderbilt)  %>% data.frame()
  
  return(Res)
}


Metadata = function(x){
  x=t(x)
  Response = rep(NA, nrow(x)) %>%  magrittr::set_names(rownames(x))
  x$patient = str_extract(rownames(x), pattern = '^ART\\d{2,3}')
  Response[which(x$patient %in% c('ART44','ART45','ART46','ART100'))] ='RCBIII'
  Response[which(x$patient %in% c('ART02','ART13','ART21','ART43','ART55','ART63','ART75','ART77','ART89'))] ='RCBII'
  Response[which(x$patient %in% c('ART10','ART17','ART51','ART66','ART70','ART92','ART94','ART134'))] = 'RCBI'
  Response[which(x$patient %in% c('ART05','ART07','ART12','ART18','ART19','ART23','ART30','ART31','ART50','ART53','ART60','ART62',
                                  'ART71','ART90','ART96','ART104','ART124'))] = 'PCR'
  Response[is.na(Response)] = 'Unknown'
  print(table(Response))
  
  GenResponse = rep(NA, length(x$patient)) %>% magrittr::set_names(colnames(x))
  GenResponse[which(x$patient %in% c('ART44','ART45','ART46','ART100',
                                     'ART02','ART13','ART21','ART43','ART55','ART63','ART75','ART77','ART89', 
                                     'ART10','ART17','ART51','ART66','ART70','ART92','ART94','ART134'))] ='NonResponder'
  GenResponse[which(x$patient %in% c('ART05','ART07','ART12','ART18','ART19','ART23','ART30','ART31','ART50','ART53','ART60','ART62',
                                     'ART71','ART90','ART96','ART104','ART124'))] = 'Responder'
  GenResponse[is.na(GenResponse)] = 'Unknown'
  print(table(GenResponse))
  
  Chem = rep(NA, length(x$patient)) %>% magrittr::set_names(colnames(x))
  Chem[which(x$patient %in% c('ART02','ART05','ART07','ART10','ART12','ART13','ART18','ART19','ART44','ART62'))] ='V2'
  Chem[is.na(Chem)] = 'V3'
  print(table(Chem))
  
  Vanderbilt = rep(NA, length(x$patient)) %>% magrittr::set_names(colnames(x))
  Vanderbilt[which(x$patient %in% c('ART17','ART18','ART44','ART43','ART50','ART70'))] ='BL1'
  Vanderbilt[which(x$patient %in% c('ART08','ART62','ART71','ART100','ART104'))] ='BL2'
  Vanderbilt[which(x$patient %in% c('ART02','ART07','ART21','ART19','ART23','ART30','ART31',
                                    'ART39','ART53','ART60','ART78','ART90','ART125','ART156','ART163','ART179','ART181'))] ='IM'
  Vanderbilt[which(x$patient %in% c('ART89','ART174'))] ='LAR'
  Vanderbilt[which(x$patient %in% c('ART05','ART12','ART13','ART46','ART51','ART55','ART66','ART75','ART92','ART94','ART124','ART134','ART170'))] ='M'
  Vanderbilt[which(x$patient %in% c('ART45'))] ='MSL'
  Vanderbilt[is.na(Vanderbilt)] = 'UNS'
  print(table(Vanderbilt))
  
  Res=cbind(Response,GenResponse,Chem,Vanderbilt)  %>% data.frame()
  
  return(Res)
}

PCA = function(mat,nPC){
  var_mad.pca <- prcomp(t(mat),
                        center = TRUE,
                        scale. = TRUE) 
  var_mad.pca.50 = var_mad.pca$x[,1:nPC]
  return(var_mad.pca.50)
}


mycomputeCommunProb = function (object, LR.use = NULL, trim = NULL, nboot = 100, seed.use = 1L, 
                                Kh = 0.5) 
{
  data <- object@data.project
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  else {
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  my.sapply <- ifelse(test = future::nbrOfWorkers() == 1, yes = pbapply::pbsapply, 
                      no = future.apply::future_sapply)
  ptm = Sys.time()
  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- length(unique(group))
  data.use <- data/max(data)
  nC <- ncol(data.use)
  index.singleL <- which(geneL %in% rownames(data.use))
  dataL1 <- data.use[geneL[index.singleL], ]
  dataL <- matrix(nrow = nLR, ncol = nC)
  dataL[index.singleL, ] <- dataL1
  index.complexL <- setdiff(1:nLR, index.singleL)
  if (length(index.complexL) > 0) {
    complex <- geneL[index.complexL]
    data.complex <- computeExpr_complex(complex_input, data.use, 
                                        complex)
    dataL[index.complexL, ] <- data.complex
  }
  index.singleR <- which(geneR %in% rownames(data.use))
  dataR1 <- data.use[geneR[index.singleR], ]
  dataR <- matrix(nrow = nLR, ncol = nC)
  dataR[index.singleR, ] <- dataR1
  index.complexR <- setdiff(1:nLR, index.singleR)
  if (length(index.complexR) > 0) {
    complex <- geneR[index.complexR]
    data.complex <- computeExpr_complex(complex_input, data.use, 
                                        complex)
    dataR[index.complexR, ] <- data.complex
  }
  index.co.A.receptor <- which(!is.na(pairLRsig$co_A_receptor) & 
                                 pairLRsig$co_A_receptor != "")
  if (length(index.co.A.receptor) > 0) {
    dataR.co.A.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                  data.use, coreceptor.all = pairLRsig$co_A_receptor, 
                                                  index.coreceptor = index.co.A.receptor)
  }
  else {
    dataR.co.A.receptor <- matrix(1, nrow = nrow(dataR), 
                                  ncol = ncol(dataR))
  }
  index.co.I.receptor <- which(!is.na(pairLRsig$co_I_receptor) & 
                                 pairLRsig$co_I_receptor != "")
  if (length(index.co.I.receptor) > 0) {
    dataR.co.I.receptor <- computeExpr_coreceptor(cofactor_input, 
                                                  data.use, coreceptor.all = pairLRsig$co_I_receptor, 
                                                  index.coreceptor = index.co.I.receptor)
  }
  else {
    dataR.co.I.receptor <- matrix(1, nrow = nrow(dataR), 
                                  ncol = ncol(dataR))
  }
  dataR <- dataR * dataR.co.A.receptor/dataR.co.I.receptor
  if (is.null(trim)) {
    FunMean <- triMean
  }
  else {
    FunMean <- function(x) mean(x, trim = trim, na.rm = TRUE)
  }
  dataLavg <- aggregate(t(dataL), list(group), FUN = FunMean)
  dataLavg <- t(dataLavg[, -1])
  dataRavg <- aggregate(t(dataR), list(group), FUN = FunMean)
  dataRavg <- t(dataRavg[, -1])
  dataL.binary = (dataL > 0) * 1
  dataR.binary = (dataR > 0) * 1
  dataLavg2 <- aggregate(t(dataL.binary), list(group), FUN = sum)
  dataLavg2 <- t(dataLavg2[, -1])/nC
  dataRavg2 <- aggregate(t(dataR.binary), list(group), FUN = sum)
  dataRavg2 <- t(dataRavg2[, -1])/nC
  index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != 
                           "")
  index.antagonist <- which(!is.na(pairLRsig$antagonist) & 
                              pairLRsig$antagonist != "")
  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  Prob <- array(0, dim = c(numCluster, numCluster, nLR))
  Pval <- array(0, dim = c(numCluster, numCluster, nLR))
  for (i in 1:nLR) {
    dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], nrow = 1), 
                                matrix(dataRavg[i, ], nrow = 1))
    P1 <- dataLR/(Kh + dataLR)
    if (sum(P1) == 0) {
      Pnull = P1
      Prob[, , i] <- Pnull
      p = 1
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                            byrow = FALSE)
    }
    else {
      if (is.element(i, index.agonist)) {
        data.agonist <- computeExpr_agonist(data.use = data.use, 
                                            pairLRsig, cofactor_input, group = group, index.agonist = i, 
                                            Kh = Kh, FunMean = FunMean)
        P2 <- Matrix::crossprod(matrix(data.agonist, 
                                       nrow = 1))
      }
      else {
        P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      if (is.element(i, index.antagonist)) {
        data.antagonist <- computeExpr_antagonist(data.use = data.use, 
                                                  pairLRsig, cofactor_input, group = group, index.antagonist = i, 
                                                  Kh = Kh, FunMean = FunMean)
        P3 <- Matrix::crossprod(matrix(data.antagonist, 
                                       nrow = 1))
      }
      else {
        P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      P4 <- Matrix::crossprod(matrix(dataLavg2[i, ], nrow = 1), 
                              matrix(dataRavg2[i, ], nrow = 1))
      Pnull = P1 * P2 * P3 * P4
      Prob[, , i] <- Pnull
      Pnull <- as.vector(Pnull)
      dataL.i <- dataL[i, ]
      dataR.i <- dataR[i, ]
      dataL2.i <- dataL.binary[i, ]
      dataR2.i <- dataR.binary[i, ]
      Pboot <- my.sapply(X = 1:nboot, FUN = function(nE) {
        groupboot <- group[permutation[, nE]]
        dataLavgB <- aggregate(matrix(dataL.i, ncol = 1), 
                               list(groupboot), FUN = FunMean)
        dataLavgB <- t(dataLavgB[, -1])
        dataLavgB <- matrix(dataLavgB, nrow = 1)
        dataRavgB <- aggregate(matrix(dataR.i, ncol = 1), 
                               list(groupboot), FUN = FunMean)
        dataRavgB <- t(dataRavgB[, -1])
        dataRavgB <- matrix(dataRavgB, nrow = 1)
        dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
        P1.boot <- dataLRB/(Kh + dataLRB)
        if (is.element(i, index.agonist)) {
          data.agonist <- computeExpr_agonist(data.use = data.use, 
                                              pairLRsig, cofactor_input, group = groupboot, 
                                              index.agonist = i, Kh = Kh, FunMean = FunMean)
          P2.boot <- Matrix::crossprod(matrix(data.agonist, 
                                              nrow = 1))
        }
        else {
          P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        if (is.element(i, index.antagonist)) {
          data.antagonist <- computeExpr_antagonist(data.use = data.use, 
                                                    pairLRsig, cofactor_input, group = groupboot, 
                                                    index.antagonist = i, Kh = Kh, FunMean = FunMean)
          P3.boot <- Matrix::crossprod(matrix(data.antagonist, 
                                              nrow = 1))
        }
        else {
          P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
        }
        dataLavg2B <- by(matrix(dataL2.i, ncol = 1), 
                         groupboot, sum)/nC
        dataLavg2B <- matrix(dataLavg2B, nrow = 1)
        dataRavg2B <- by(matrix(dataR2.i, ncol = 1), 
                         groupboot, sum)/nC
        dataRavg2B <- matrix(dataRavg2B, nrow = 1)
        P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
        Pboot = P1.boot * P2.boot * P3.boot * P4.boot
        return(as.vector(Pboot))
      })
      Pboot <- matrix(unlist(Pboot), nrow = length(Pnull), 
                      ncol = nboot, byrow = FALSE)
      nReject <- rowSums(Pboot - Pnull >= 0)
      p = nReject/nboot
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, 
                            byrow = FALSE)
    }
  }
  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list(prob = Prob, pval = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")
  object@net <- net
  return(object)
}
environment(mycomputeCommunProb) <- environment(computeCommunProb)
plotPCA.san <- function (object, intgroup = "condition", PC1=1,PC2=2, ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC1 = pca$x[, PC1], PC2 = pca$x[, PC2], group = group, 
                  intgroup.df, name = colData(rld)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(PC1,PC2)]
    return(d)
  }
  ggplot(data = d, aes_string(x ='PC1', y = 'PC2', color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC",PC1,": ", round(percentVar[PC1] * 100), "% variance")) + ylab(paste0("PC",PC2,": ", round(percentVar[PC2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

readVarbinCNA_copykat<-function (dir, remove_Y = FALSE, genome_version = c("hg19"), 
                                 bin_size = c("200k", "100k"), gamma =40, clean_names = TRUE,method = NULL) 
{
  if (fs::file_exists(fs::dir_ls(path = dir, recurse = T, glob = "*uber*seg.txt")) == 
      FALSE) {
    stop("Segment ratio matrix can't be found in the provided directory.\n      Please make sure a uber.seg file can be found.")
  }
  if (length(fs::dir_ls(path = dir, recurse = T, glob = "*uber*seg.txt")) > 
      1) {
    stop("More than one uber.seg file can be found at the provided directory.\n      Please make sure to only have one sample at that location.")
  }
  message("Importing segment ratios.")
  dat <- data.table::fread(input = fs::dir_ls(path = dir, recurse = T, 
                                              glob = "*uber*seg.txt"), showProgress = TRUE, integer64 = "double") %>% 
    as.data.frame()
  dim(dat)
  dat$chrom %>%table()
  rg = dat[,c(1:3)]
  dat <- t(dat[, -c(1:3)]) %>% set_colnames(paste0('chr', dat$chrom, '_', dat$chrompos)) %>% t()
  dat[1:5,1:5]
  dat = cbind(rg,dat)
  
  if (clean_names == TRUE) {
    dat <- janitor::clean_names(dat)
  }
  if (remove_Y == TRUE) {
    dat <- dat %>% dplyr::filter(chrom != 24)
  }
  seg_data <- dat %>% dplyr::select(-c(chrom, chrompos, abspos))
  message("Importing ratios.")
  dat_rat <- seg_data
  message("Importing bin counts.")
  dat_bin <- seg_data
  seg_data[1:5,1:5]
  
  ref = hg19_rg[1:nrow(dat),]
  ref_chrarm <- ref %>%
    dplyr::mutate(chrarm = paste0(gsub("chr", "", chr), arm))
  levels_chrarm <- gtools::mixedsort(unique(ref_chrarm$chrarm))
  ref_chrarm <- ref_chrarm %>%
    dplyr::mutate(chrarm = as.factor(chrarm)) %>%
    dplyr::mutate(chrarm = forcats::fct_relevel(chrarm, levels_chrarm))
  ref_chrarm = ref_chrarm
  seg_df = seg_data
  
  if (method == 'multipcf'){
    mpcf <- copynumber::multipcf(dat %>% select(-abspos),
                                 arms = vapply(
                                   regmatches(ref_chrarm$chrarm,
                                              regexec("[pq]",
                                                      ref_chrarm$chrarm)),
                                   FUN =  "[",
                                   1,
                                   FUN.VALUE = character(1)
                                 )
    )
    
    seg_df <- apply(mpcf[, 6:ncol(mpcf)], 2, function(x) {
      rep.int(x, mpcf$n.probes)
    })
    seg_df <- round(as.data.frame(seg_df), 2)
    seg_df[1:5,1:5]
  }
  
  if (method == "CBS") {
    seg_list <-
      BiocParallel::bplapply(
        seg_data,
        FUN = function(x) {
          CNA_object <-
            DNAcopy::CNA(x,
                         ref_chrarm$chrarm,
                         ref$start,
                         data.type = "logratio",
                         sampleid = names(x)
            )
          
          withr::with_seed(seed = 17,
                           segment_smoothed_CNA_object <-
                             .quiet(
                               DNAcopy::segment(
                                 CNA_object,
                                 alpha = 1e-5,
                                 min.width = 5,
                                 undo.splits = 'prune'
                               )
                             )
          )
          
          
          short_cbs <- segment_smoothed_CNA_object[[2]]
          log_seg_mean_LOWESS <-
            rep(short_cbs$seg.mean, short_cbs$num.mark)
        },
        BPPARAM = bpparam()
      )
    
    seg_df <- dplyr::bind_cols(seg_list) %>%
      as.data.frame() %>%
      round(2)
  }
  
  if (method == 'origin') {
    seg_df = seg_data
  }
  rownames(seg_df) = rownames(seg_data)
  
  #set range
  if(genome_version == 'hg19'){
    gr_varbin_full <-makeGRangesFromDataFrame(hg19_rg,keep.extra.columns = T)
  }else{
    gr_varbin_full <-makeGRangesFromDataFrame(hg38_rg,keep.extra.columns = T)
    
  }
  bin_size <- '200k'
  # grlist_varbin <- switch(genome_version, hg19 = varbin_hg19_grangeslist)
  # tmp_key <- paste0("res_", bin_size)
  # gr_varbin_full <- grlist_varbin[[tmp_key]]
  GenomeInfoDb::seqlevelsStyle(gr_varbin_full) <- "Ensembl"
  gr_varbin_full <- GenomeInfoDb::renameSeqlevels(gr_varbin_full,  c(X = 23, Y = 24))
  
  rg <- dat %>% dplyr::select(c(chrom, chrompos, abspos)) %>% as.data.frame()
  IRanges::start(gr_varbin_full) %>% length()
  length(rg$chrompos)
  
  key_ref <- paste0(GenomicRanges::seqnames(gr_varbin_full), 
                    "_", IRanges::start(gr_varbin_full))
  idx = grepl('24_',x=key_ref)
  
  g <- gr_varbin_full[!idx, ]
  g$abspos <- rg$abspos
  g <- GenomeInfoDb::renameSeqlevels(g, c(`23` = "X", `24` = "Y"))
  GenomeInfoDb::seqlevelsStyle(g) <- "UCSC"
  cna_obj <- scCNA(segment_ratios = seg_df, ratios = dat_rat, 
                   bin_counts = dat_bin, rowRanges = g)
  SummarizedExperiment::colData(cna_obj)$sample <- names(seg_df)
  if (rlang::is_empty(fs::dir_ls(path = dir, recurse = T, glob = "*stat_metrics.txt"))) {
    warning("No metrics file found. \n\n            Metrics are needed if you'd like to run copykit::runMetrics()\n\n            Make sure folder metrics with file all_stat_metrics.txt can be found by copykit::runVarbinCNA()")
  }else {
    if (fs::file_exists(fs::dir_ls(path = dir, recurse = T, 
                                   glob = "*stat_metrics.txt"))) {
      message("Importing metrics.")
      dat_metrics <- data.table::fread(fs::dir_ls(path = dir, 
                                                  recurse = T, glob = "*stat_metrics.txt"), showProgress = TRUE, 
                                       integer64 = "double") %>% janitor::clean_names() %>% 
        dplyr::rename(sample = "sample_name") %>% as.data.frame()
      if (clean_names == TRUE) {
        dat_metrics <- dat_metrics %>% dplyr::mutate(sample = janitor::make_clean_names(sample))
      }
      dat_metrics <- dat_metrics[match(SummarizedExperiment::colData(cna_obj)$sample, 
                                       dat_metrics$sample), ]
      if (identical(dat_metrics$sample, SummarizedExperiment::colData(cna_obj)$sample)) {
        SummarizedExperiment::colData(cna_obj)$reads_total <- dat_metrics$total_reads
        SummarizedExperiment::colData(cna_obj)$reads_assigned_bins <- dat_metrics$reads_kept
        SummarizedExperiment::colData(cna_obj)$percentage_duplicates <- round(dat_metrics$dups_removed/dat_metrics$total_reads,2)
      }else {
        a = SummarizedExperiment::colData(cna_obj)$sample
        dat_metrics$sample = a
        SummarizedExperiment::colData(cna_obj)$reads_total <- dat_metrics$total_reads
        SummarizedExperiment::colData(cna_obj)$reads_assigned_bins <- dat_metrics$reads_kept
        SummarizedExperiment::colData(cna_obj)$percentage_duplicates <- round(dat_metrics$dups_removed/dat_metrics$total_reads,2)
      }
    }
  }
  if (remove_Y == TRUE) {
    message("Removed ChrY information.")
  }
  return(cna_obj)
}


#################### scCNA for copykit ########################3
# scCNA = function (segment_ratios, ratios, bin_counts, consensus = data.frame(), 
#           phylo = structure(list(), class = "phylo"), consensusPhylo = structure(list(), 
#                                                                                  class = "phylo"), distMat = dist(matrix(0, 0, 0)), graph = igraph::graph.empty(), 
#           ...) {
#   cna <- SingleCellExperiment::SingleCellExperiment(list(segment_ratios = segment_ratios, 
#                                                          ratios = ratios, bin_counts = bin_counts), ...)
#   .scCNA(cna, phylo = phylo, consensusPhylo = consensusPhylo, 
#          distMat = distMat, graph = graph, consensus = consensus)
# }


plotHeatmap_copykat = function (scCNA, assay = "segment_ratios", pt.name = '',
                                order_cells = c("consensus_tree", "hclust", "phylogeny"), 
                                label = NULL, label_colors = NULL, lowscale = -0.5,highscale = 0.5,
                                consensus = FALSE, rounding_error = FALSE, row_split = NULL, n_threads = 48) 
{
  order_cells <- match.arg(order_cells)
  if (is.null(label) & !is.null(label_colors)) {
    stop("Please provide a label argument if colors are being specified for it.")
  }
  if (!is.null(label_colors) & !is.list(label_colors)) {
    stop("label_colors argument must be a named list.")
  }
  if (!is.null(label_colors)) {
    if (length(label_colors) != length(label)) {
      stop("Label and Label colors arguments must have the same length.")
    }
  }
  if (!is.null(label) & !is.character(label)) {
    stop("Label must be a character vector.")
  }
  if (is.null(SummarizedExperiment::colData(scCNA)$subclones) && 
      order_cells != "phylogeny") {
    message("Ordering by consensus requires cluster information.\nSwitching to hclust")
    order_cells <- "hclust"
  }
  if (rounding_error == TRUE && assay != "integer") {
    stop("Rounding error argument must be used with assay 'integer'.")
  }
  if (consensus == TRUE) {
    if (attr(consensus(scCNA), "consensus_assay") == "integer") {
      assay <- "integer"
    }
  }
  
  color_heat = circlize::colorRamp2(breaks = c(lowscale,0,highscale), c("dodgerblue3", "white", "firebrick3"))
  if (assay == "integer") {
    mean_ploidy <- mean(SummarizedExperiment::colData(scCNA)$ploidy)
    ploidy_trunc <- 2 * round(mean_ploidy)
    color_heat <- structure(pals::ocean.balance(length(0:ploidy_trunc)), 
                            names = 0:ploidy_trunc)
    if (round(mean_ploidy) == 2) {
      color_heat <- structure(c("#3787BA", "#95B8C5", "#F0ECEB", 
                                "#D7A290", "#BF583B", "#8D1128", "#3C0912"), 
                              names = c("0", "1", "2", "3", "4", "5", "6"))
    }
  }
  
  color_heat = circlize::colorRamp2(breaks = c(lowscale,0,highscale), c("dodgerblue3", "white", "firebrick3"))
  
  seg_data <- t(SummarizedExperiment::assay(scCNA, assay))
  seg_data[1:5,1:5]
  chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
  if (any(chr_ranges$seqnames == "24") || any(chr_ranges$seqnames == 
                                              "Y") || any(chr_ranges$seqnames == "chrY")) {
    chr_binary <- rep(c(2, 1), length(chr_lengths)/2)
  }
  else {
    chr_binary <- c(rep(c(2, 1), (length(chr_lengths)/2)), 
                    2)
  }
  chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  chr_rl_c <- c(1, cumsum(chr_lengths))
  chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], 
                       b = chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))
  chrom.names <- c(1:22, "X", "Y")
  v <- vector(length = sum(chr_lengths), mode = "character")
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""
  chr_bar <- ComplexHeatmap::HeatmapAnnotation(chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)], 
                                                                                    gp = grid::gpar(fontsize = 14)), df = as.character(chr[1:nrow(chr), 
                                                                                    ]), show_legend = FALSE, show_annotation_name = FALSE, 
                                               which = "column", col = list(df = c(`1` = "grey88", `2` = "black")))
  if (consensus == FALSE) {
    if (order_cells == "consensus_tree") {
      if (nrow(consensus(scCNA)) == 0) {
        scCNA <- calcConsensus(scCNA)
      }
      scCNA <- runConsensusPhylo(scCNA)
      consensus_by <- attr(consensus(scCNA), "consensus_by")
      meta <- as.data.frame(colData(scCNA)) %>% dplyr::select(sample, 
                                                              !!consensus_by)
      meta_info <- as.character(dplyr::pull(meta, !!consensus_by))
      tree <- consensusPhylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      meta_o <- meta[order(match(meta_info, tree_tips_order)), 
      ]
      seg_data_ordered <- seg_data[meta_o$sample, ]
    }else{
      seg_data_ordered <- seg_data
    }
    if (order_cells == "phylogeny") {
      if (ape::Ntip(phylo(scCNA)) == 0) {
        stop("No phylogeny detected in scCNA object. Use runPhylo")
      }
      tree <- phylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      seg_data_ordered <- seg_data[tree_tips_order, ]
    }
    if (order_cells == "hclust") {
      if (length(copykit::distMat(scCNA)) == 0) {
        message("No distance matrix detected in the scCNA object.")
        scCNA <- runDistMat(scCNA, metric = "euclidean",n_threads = n_threads)
      }
      if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
        stop("Number of samples in the distance matrix different\n          from number of samples in the scCNA object.\n        Perhaps you filtered your dataset? use copykit::runDistMat() to update it.")
      }
      hc <- fastcluster::hclust(distMat(scCNA), method = "ward.D2")
      seg_data_ordered <- seg_data[hc$order, ]
    }
  }
  else {
    scCNA <- runConsensusPhylo(scCNA)
    tree <- consensusPhylo(scCNA)
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
      rev()
    seg_data <- as.matrix(t(consensus(scCNA)))
    seg_data_ordered <- seg_data[tree_tips_order, ]
  }
  if (assay == "integer") {
    seg_data_int <- seg_data_ordered
    seg_data_int[seg_data_int > ploidy_trunc] <- ploidy_trunc
    if (rounding_error == TRUE) {
      int_nr <- as.matrix(SummarizedExperiment::assay(scCNA, 
                                                      "segment_ratios")) %*% diag(SummarizedExperiment::colData(scCNA)$ploidy)
      names(int_nr) <- names(SummarizedExperiment::assay(scCNA, 
                                                         "segment_ratios"))
      err <- abs(int_nr - assay(scCNA, "integer"))
      err <- err[rownames(seg_data_int)]
    }
  }
  message("Plotting Heatmap.")
  complex_args <- list(use_raster = TRUE, column_title = "genomic coordinates", 
                       column_title_gp = grid::gpar(fontsize = 18), column_title_side = "bottom", 
                       row_title = paste0(pt.name, '             ',nrow(seg_data_ordered), " samples"), 
                       row_title_gp = grid::gpar(fontsize = 18), top_annotation = chr_bar, 
                       cluster_rows = FALSE, border = TRUE, cluster_columns = FALSE, 
                       show_column_names = FALSE, show_row_names = FALSE, show_heatmap_legend = TRUE)
  if (is.null(label)) {
    if (assay != "integer") {
      suppressMessages(do.call(ComplexHeatmap::Heatmap, 
                               c(list(matrix = seg_data_ordered, 
                                      heatmap_legend_param = list(title = "log2 (segratio)"),
                                      col = color_heat), 
                                 complex_args)))
    }
    else {
      if (rounding_error == FALSE) {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                heatmap_legend_param = list(title = "copy number"), 
                                                col = color_heat), complex_args))
      }
      else {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                heatmap_legend_param = list(title = "rounding error"), 
                                                col = viridis::viridis(200)), complex_args))
      }
    }
  }
  else {
    metadata <- SummarizedExperiment::colData(scCNA) %>% 
      as.data.frame()
    if (consensus == FALSE) {
      metadata <- metadata[rownames(seg_data_ordered), 
      ]
    }
    metadata_anno_df <- metadata %>% dplyr::select(dplyr::all_of(label))
    if (consensus == TRUE) {
      cons_attr <- attr(consensus(scCNA), "consensus_by")
      if (length(label) > 1) {
        stop("Label must be of length 1 for consensus heatmap annotation.")
      }
      if (cons_attr != label) {
        stop("Consensus heatmap can only be annotated with the same metadata element\n             used for generating the consensus matrix.")
      }
      metadata_anno_df <- metadata_anno_df[label] %>% dplyr::distinct()
      rownames(metadata_anno_df) <- metadata_anno_df %>% 
        dplyr::pull(!!cons_attr)
      metadata_anno_df <- metadata_anno_df[rownames(seg_data_ordered), 
                                           , drop = FALSE]
    }
    if (is.null(label_colors)) {
      h <- 15
      l <- 65
      cont_options <- c("D", "A", "E", "C", "B")
      cont_i <- 1
      label_colors <- list()
      label_colors <- c(list(superclones = superclones_pal(), 
                             subclones = subclones_pal(), filtered = c(removed = "#DA614D", 
                                                                       kept = "#5F917A"), is_normal = c(`TRUE` = "#396DB3", 
                                                                                                        `FALSE` = "#11181D")), label_colors)
      default_labels <- c("superclones", "subclones", "filtered", 
                          "is_normal")
      for (i in 1:length(label)) {
        if (any(str_detect(label[i], default_labels))) {
          label_colors[i] <- label_colors[default_labels[stringr::str_detect(label[i], 
                                                                             default_labels)]]
          names(label_colors)[i] <- label[i]
        }
        else if (is.numeric(dplyr::pull(metadata_anno_df, 
                                        label[i]))) {
          n = 300
          min_v = min(dplyr::pull(metadata_anno_df, label[i]))
          max_v = max(dplyr::pull(metadata_anno_df, label[i]))
          label_colors[i] <- list(circlize::colorRamp2(seq(min_v, 
                                                           max_v, length = n), viridis::viridis(n, option = cont_options[cont_i])))
          names(label_colors)[i] <- label[i]
          cont_i <- cont_i + 1
        }
        else {
          elements <- metadata_anno_df %>% dplyr::pull(label[i]) %>% 
            unique() %>% as.character() %>% sort()
          n <- length(elements)
          hex <- (scales::hue_pal(h = c(0, 360) + h, 
                                  l = 65))(n)
          col <- structure(hex, names = elements)
          label_colors[i] <- list(col)
          names(label_colors)[i] <- label[i]
          l <- l - 10
          h <- h + 15
        }
      }
    }
    label_colors[sapply(label_colors, is.null)] <- NULL
    cluster_anno <- ComplexHeatmap::rowAnnotation(df = metadata_anno_df, 
                                                  col = label_colors, show_annotation_name = T)
    if (!is.null(row_split)) {
      if (length(row_split) > 1) {
        stop("row_split length must be 1")
      }
      else {
        if (assay != "integer") {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_ordered, row_split = dplyr::pull(metadata_anno_df, 
                                                                                                     row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "log2 (segratio)"),col = color_heat), 
                                             complex_args))
        }
        else {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  row_split = dplyr::pull(metadata_anno_df, 
                                                                          row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
      }
    }
    else {
      if (assay != "integer") {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_ordered, left_annotation = cluster_anno, heatmap_legend_param = list(title = "log2 (segratio)"),
                                                col = color_heat), 
                                           complex_args))
      }
      else {
        if (rounding_error == FALSE) {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
        else {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "rounding error"), 
                                                  col = viridis::viridis(200)), complex_args))
        }
      }
    }
  }
}



#################### MarkovHC_yl ########################3
MarkovHC_yl = function (MarkovHC_input = NULL, SNNslot = NULL, KNNslot = NULL, ncore =50,
                        
                        KNN = 20, dobasecluster = TRUE, cutpoint = 0.001, verbose = TRUE) 
{
  if ((!is.matrix(MarkovHC_input)) & (!(class(MarkovHC_input) == 
                                        "Seurat"))) {
    print("The type of 'MarkovHC_input' must be matirx or seurat object!")
    return(NULL)
  }
  if (is.matrix(MarkovHC_input)) {
    print("The input is a matrix.")
    origin_matrix <- t(MarkovHC_input)
  }
  else {
    print("The input is a Seurat object.")
  }
  if (!is.numeric(KNN)) {
    print("The type of 'KNN' must be numeric!")
    return(NULL)
  }
  if (!is.numeric(cutpoint)) {
    print("The type of 'cutpoint' must be numeric!")
    return(NULL)
  }
  set.seed(1)
  cl <- makeCluster(getOption("cl.cores", ncore))
  registerDoParallel(cl)
  if (is.matrix(MarkovHC_input)) {
    symmetric_KNN_graph <- FindNeighbors(object = origin_matrix, 
                                         k.param = KNN, compute.SNN = TRUE, prune.SNN = 0)
    symmetric_KNN_graph <- Seurat::as.sparse(symmetric_KNN_graph$snn) * 
      Seurat::as.sparse(symmetric_KNN_graph$nn) * Matrix::t(Seurat::as.sparse(symmetric_KNN_graph$nn))
  }
  else {
    symmetric_KNN_graph <- Seurat::as.sparse(MarkovHC_input@graphs[[SNNslot]]) * 
      Seurat::as.sparse(MarkovHC_input@graphs[[KNNslot]]) * 
      Matrix::t(Seurat::as.sparse(MarkovHC_input@graphs[[KNNslot]]))
  }
  gc(verbose = FALSE)
  rownames(symmetric_KNN_graph) <- as.character(1:nrow(symmetric_KNN_graph))
  colnames(symmetric_KNN_graph) <- as.character(1:ncol(symmetric_KNN_graph))
  symmetric_KNN_graph_object <- graph_from_adjacency_matrix(adjmatrix = symmetric_KNN_graph, 
                                                            mode = "undirected", weighted = TRUE, diag = TRUE)
  centrality_scores <- igraph::degree(symmetric_KNN_graph_object, 
                                      v = V(symmetric_KNN_graph_object), mode = "total", loops = TRUE, 
                                      normalized = FALSE)
  if (dobasecluster == TRUE) {
    cluster_louvain_object <- cluster_louvain(graph = symmetric_KNN_graph_object)
    hresult_cut <- cluster_louvain_object$memberships[1, 
    ]
  }
  else {
    hresult_cut <- 1:nrow(symmetric_KNN_graph)
  }
  unique_clusters <- unique(hresult_cut)
  if (dobasecluster == TRUE) {
    symmetric_KNN_graph_merged <- as(Matrix::Matrix(data = 0, 
                                                    nrow = length(unique_clusters), ncol = ncol(symmetric_KNN_graph), 
                                                    sparse = TRUE), "dgCMatrix")
    for (clusterindex in unique_clusters) {
      temp_index <- which(hresult_cut == clusterindex)
      if (length(temp_index) == 1) {
        symmetric_KNN_graph_merged[clusterindex, ] <- symmetric_KNN_graph[temp_index, 
        ]
      }
      else {
        temp_cluster <- symmetric_KNN_graph[temp_index, 
        ]
        symmetric_KNN_graph_merged[clusterindex, ] <- as.vector(qlcMatrix::colMax(temp_cluster, 
                                                                                  which = FALSE, ignore.zero = FALSE))
      }
    }
    symmetric_KNN_graph_cluster <- matrix(0, nrow(symmetric_KNN_graph_merged), 
                                          nrow(symmetric_KNN_graph_merged))
    rownames(symmetric_KNN_graph_cluster) <- as.character(1:nrow(symmetric_KNN_graph_cluster))
    colnames(symmetric_KNN_graph_cluster) <- as.character(1:ncol(symmetric_KNN_graph_cluster))
    for (clusterindex in unique_clusters) {
      for (clusterindex2 in unique_clusters) {
        temp_cluster <- symmetric_KNN_graph_merged[clusterindex, 
                                                   which(hresult_cut == clusterindex2)]
        symmetric_KNN_graph_cluster[clusterindex, clusterindex2] <- max(temp_cluster)
      }
    }
    diag(symmetric_KNN_graph_cluster) <- 1
    centrality_scores_cluster <- integer(length = length(unique_clusters))
    for (score_index in unique_clusters) {
      centrality_scores_tpm <- centrality_scores[which(hresult_cut == 
                                                         score_index)]
      centrality_scores_cluster[score_index] <- mean(centrality_scores_tpm)
    }
  }
  else {
    symmetric_KNN_graph_cluster <- as.matrix(symmetric_KNN_graph)
    rownames(symmetric_KNN_graph_cluster) <- as.character(1:nrow(symmetric_KNN_graph_cluster))
    colnames(symmetric_KNN_graph_cluster) <- as.character(1:ncol(symmetric_KNN_graph_cluster))
    centrality_scores_cluster <- centrality_scores
  }
  D.matrix <- matrix(0, nrow = ncol(symmetric_KNN_graph_cluster), 
                     ncol = ncol(symmetric_KNN_graph_cluster))
  rownames(D.matrix) <- rownames(symmetric_KNN_graph_cluster)
  colnames(D.matrix) <- rownames(symmetric_KNN_graph_cluster)
  nothing <- foreach(nrow_i = 1:nrow(D.matrix), .combine = "c") %do% 
    {
      D.matrix[nrow_i, ] <- (centrality_scores_cluster[nrow_i]/centrality_scores_cluster)^2
      NULL
    }
  transitionMatrix_pseudoenergyMatrix <- transition_probability_pseudo_energy_matrix(matrix = symmetric_KNN_graph_cluster, 
                                                                                     D.matrix = D.matrix)
  transitionMatrix <- transitionMatrix_pseudoenergyMatrix$P
  C_matrix <- transitionMatrix_pseudoenergyMatrix$C
  C_matrix_temp <- C_matrix
  C_matrix_temp[which(abs(C_matrix_temp) < exp(-30))] <- Inf
  ground_energy <- min(C_matrix_temp)/100
  C_matrix <- C_matrix + ground_energy
  C_matrix[which(is.infinite(C_matrix) == TRUE)] <- 0
  rownames(C_matrix) <- as.character(1:nrow(C_matrix))
  colnames(C_matrix) <- as.character(1:ncol(C_matrix))
  C_matrix_graph_object <- graph_from_adjacency_matrix(adjmatrix = as.matrix(C_matrix), 
                                                       mode = "directed", weighted = TRUE, diag = TRUE)
  if (verbose == TRUE) {
    print("Calculate the shortest distance between each vertex pair in the graph.")
  }
  C_matrix_graph_shortest_distance <- igraph::distances(graph = C_matrix_graph_object, 
                                                        v = 1:nrow(C_matrix), to = 1:nrow(C_matrix), mode = "out", 
                                                        weights = igraph::E(C_matrix_graph_object)$weight, algorithm = "dijkstra")
  P_updated <- transitionMatrix
  MarkovHC_result <- list()
  attractors <- list()
  basins <- list()
  graphvertex_attractors <- list()
  graphvertex_basins <- list()
  attractorPoints <- list()
  basinPoints <- list()
  for (i in 1:length(unique_clusters)) {
    attractors <- c(attractors, list(i))
    basins <- c(basins, list(i))
    graphvertex_attractors <- c(graphvertex_attractors, list(i))
    graphvertex_basins <- c(graphvertex_basins, list(i))
    clusterPoints <- which(hresult_cut == i)
    attractorPoints <- c(attractorPoints, list(clusterPoints))
    basinPoints <- c(basinPoints, list(clusterPoints))
  }
  basinNum <- length(attractors)
  level_result <- list(basins = basins, attractors = attractors, 
                       graphvertex_attractors = graphvertex_attractors, graphvertex_basins = graphvertex_basins, 
                       basinPoints = basinPoints, attractorPoints = attractorPoints, 
                       basinNum = basinNum, energyCutpoint = 0, C_matrix_updatedmatrix = C_matrix_graph_shortest_distance)
  MarkovHC_result <- append(MarkovHC_result, list(level_result))
  levels_indice <- 1
  if (verbose == TRUE) {
    print(paste("Build the level ", as.character(levels_indice), 
                "...", sep = ""))
  }
  while (TRUE) {
    levels_indice <- levels_indice + 1
    if (verbose == TRUE) {
      print(paste("Build the level ", as.character(levels_indice), 
                  "...", sep = ""))
    }
    RS_vector <- judge_RS(P = P_updated)
    attractors <- list()
    basins <- list()
    graphvertex_attractors <- list()
    graphvertex_basins <- list()
    attractorPoints <- list()
    basinPoints <- list()
    processed_attractors <- integer(length = length(RS_vector))
    processed_attractors[which(RS_vector <= exp(-10))] <- 1
    attractor_indice <- 1
    P_updated_graph_object <- graph_from_adjacency_matrix(adjmatrix = as.matrix(P_updated), 
                                                          mode = "directed", weighted = TRUE, diag = TRUE)
    while (TRUE) {
      if (all(processed_attractors == 1)) {
        break
      }
      if (verbose == TRUE) {
        print(paste("Find attractors in the basin ", 
                    as.character(attractor_indice), ".", sep = ""))
      }
      attractor_temp <- which(processed_attractors == 0)[1]
      attractor_temp_access <- subcomponent(graph = P_updated_graph_object, 
                                            v = attractor_temp, mode = "out") %>% as.vector()
      processed_attractors[attractor_temp_access] <- 1
      attractors <- c(attractors, list(attractor_temp_access))
      basins_temp_merged <- subcomponent(graph = P_updated_graph_object, 
                                         v = attractor_temp_access, mode = "in") %>% as.vector()
      processed_attractors[basins_temp_merged] <- 1
      basins <- c(basins, list(basins_temp_merged))
      attractor_indice <- attractor_indice + 1
    }
    if (length(unique(unlist(basins))) != nrow(P_updated)) {
      if (verbose == TRUE) {
        print("Correct the lost basins caused by R precision limitation.")
      }
      lost_basin <- setdiff(1:nrow(P_updated), unique(unlist(basins)))
      lost_attractor <- integer(length = length(lost_basin))
      for (i in 1:length(lost_basin)) {
        Si <- subcomponent(graph = P_updated_graph_object, 
                           v = lost_basin[i], mode = "out") %>% as.vector()
        Ti <- subcomponent(graph = P_updated_graph_object, 
                           v = lost_basin[i], mode = "in") %>% as.vector()
        if (length(setdiff(Si, Ti)) > 0) {
          lost_attractor[i] <- 1
        }
        else {
          lost_attractor[i] <- 0
        }
      }
      while (TRUE) {
        if (all(lost_attractor == 1)) {
          break
        }
        if (verbose == TRUE) {
          print(paste("Find attractors in the basin ", 
                      as.character(attractor_indice), ".", sep = ""))
        }
        attractor_temp <- lost_basin[which(lost_attractor == 
                                             0)[1]]
        attractor_temp_access <- subcomponent(graph = P_updated_graph_object, 
                                              v = attractor_temp, mode = "out") %>% as.vector()
        lost_attractor[which(lost_basin %in% attractor_temp_access)] <- 1
        attractors <- c(attractors, list(attractor_temp_access))
        basins_temp_merged <- subcomponent(graph = P_updated_graph_object, 
                                           v = attractor_temp_access, mode = "in") %>% 
          as.vector()
        lost_attractor[which(lost_basin %in% basins_temp_merged)] <- 1
        basins <- c(basins, list(basins_temp_merged))
        attractor_indice <- attractor_indice + 1
      }
    }
    basinNum <- length(basins)
    basin_indice <- 1
    for (i in 1:basinNum) {
      if (verbose == TRUE) {
        print(paste("Partition the basin ", as.character(basin_indice), 
                    ".", sep = ""))
      }
      attractorPoints_temp <- level_result$attractorPoints[attractors[[i]]] %>% 
        unlist() %>% unique()
      attractorPoints <- c(attractorPoints, list(attractorPoints_temp))
      basinPoints_temp <- level_result$basinPoints[basins[[i]]] %>% 
        unlist() %>% unique()
      basinPoints <- c(basinPoints, list(basinPoints_temp))
      graphvertex_attractors_temp <- level_result$graphvertex_attractors[attractors[[i]]] %>% 
        unlist() %>% unique()
      graphvertex_attractors <- c(graphvertex_attractors, 
                                  list(graphvertex_attractors_temp))
      graphvertex_basins_temp <- level_result$graphvertex_basins[basins[[i]]] %>% 
        unlist() %>% unique()
      graphvertex_basins <- c(graphvertex_basins, list(graphvertex_basins_temp))
      basin_indice <- basin_indice + 1
    }
    if (verbose == TRUE) {
      print("Update the pseudo energy matrix.")
    }
    changed_basins <- which((sapply(basins, length) > 1) == 
                              TRUE)
    unchanged_basins <- which((sapply(basins, length) == 
                                 1) == TRUE)
    C_matrix_updated <- matrix(data = 0, nrow = basinNum, 
                               ncol = basinNum)
    C_matrix_updated[unchanged_basins, unchanged_basins] <- level_result$C_matrix_updatedmatrix[unlist(basins[unchanged_basins]), 
                                                                                                unlist(basins[unchanged_basins])]
    for (i in changed_basins) {
      for (j in changed_basins) {
        source_v <- level_result$graphvertex_attractors[attractors[[i]]] %>% 
          unlist() %>% unique()
        target_v <- level_result$graphvertex_basins[basins[[j]]] %>% 
          unlist() %>% unique()
        C_matrix_updated[i, j] <- C_matrix_graph_shortest_distance[source_v, 
                                                                   target_v] %>% min()
      }
    }
    for (i in changed_basins) {
      for (j in unchanged_basins) {
        source_v <- level_result$graphvertex_attractors[attractors[[i]]] %>% 
          unlist() %>% unique()
        target_v <- level_result$graphvertex_basins[basins[[j]]] %>% 
          unlist() %>% unique()
        C_matrix_updated[i, j] <- C_matrix_graph_shortest_distance[source_v, 
                                                                   target_v] %>% min()
        source_v <- level_result$graphvertex_attractors[attractors[[j]]] %>% 
          unlist() %>% unique()
        target_v <- level_result$graphvertex_basins[basins[[i]]] %>% 
          unlist() %>% unique()
        C_matrix_updated[j, i] <- C_matrix_graph_shortest_distance[source_v, 
                                                                   target_v] %>% min()
      }
    }
    C_matrix_updated_indice <- C_matrix_updated
    diag(C_matrix_updated_indice) <- Inf
    if ((basinNum > 1) & (all(is.infinite(C_matrix_updated_indice)) == 
                          FALSE)) {
      if (verbose == TRUE) {
        print("Update the transition probability matrix.")
      }
      update_P_result <- update_P(C_matrix_updated = C_matrix_updated, 
                                  C_cut = cutpoint)
      P_updated <- update_P_result[[1]]
      energyCutpoint <- update_P_result[[2]]
    }
    level_result <- list(basins = basins, attractors = attractors, 
                         graphvertex_attractors = graphvertex_attractors, 
                         graphvertex_basins = graphvertex_basins, basinPoints = basinPoints, 
                         attractorPoints = attractorPoints, basinNum = basinNum, 
                         C_matrix_updatedmatrix = C_matrix_updated, P_updated = P_updated, 
                         energyCutpoint = energyCutpoint, changed_basins = changed_basins)
    MarkovHC_result <- append(MarkovHC_result, list(level_result))
    if ((basinNum == 1) | (all(is.infinite(C_matrix_updated_indice)))) {
      inputParameters <- list(KNN = KNN, dobasecluster = dobasecluster, 
                              cutpoint = cutpoint)
      C_cut_seq <- ground_energy
      for (i in 2:length(MarkovHC_result)) {
        C_cut_seq <- c(C_cut_seq, MarkovHC_result[[i]][["energyCutpoint"]])
      }
      midResults <- list(symmetric_KNN_graph = symmetric_KNN_graph, 
                         symmetric_KNN_graph_object = symmetric_KNN_graph_object, 
                         centrality_scores = centrality_scores, centrality_scores_cluster = centrality_scores_cluster, 
                         symmetric_KNN_graph_cluster = symmetric_KNN_graph_cluster, 
                         transitionMatrix = transitionMatrix, C_matrix = C_matrix, 
                         C_matrix_graph_object = C_matrix_graph_object, 
                         C_matrix_graph_shortest_distance = C_matrix_graph_shortest_distance, 
                         C_cut_seq = C_cut_seq, ground_energy = ground_energy)
      MarkovHC_object <- list(hierarchicalStructure = MarkovHC_result, 
                              inputParameters = inputParameters, midResults = midResults)
      names(MarkovHC_object$hierarchicalStructure) <- paste(rep("level", 
                                                                length(MarkovHC_result)), 1:length(MarkovHC_result), 
                                                            sep = "")
      stopCluster(cl)
      return(MarkovHC_object)
    }
  }
}

select_pc = function(SeuratObject){
  model <- SiZer::piecewise.linear(x = 1:length(SeuratObject@reductions$pca@stdev), 
                                   y = SeuratObject@reductions$pca@stdev)
  recommend_level <- ceiling((model$change.point) + 1)
  plot(model)
  return(recommend_level)
}


clust_gene_ovlp <- function(top_de, marker_list, bk_gene) {
  top_de$cluster %<>% as.character()
  clust_query <- unique(top_de$cluster)
  clust_pval_df <- matrix(NA, length(marker_list), length(clust_query)) %>%
    set_rownames(names(marker_list)) %>% set_colnames(clust_query) %>% data.frame
  clust_or_df <- matrix(NA, length(marker_list), length(clust_query)) %>%
    set_rownames(names(marker_list)) %>% set_colnames(clust_query) %>% data.frame
  for (i in 1:length(clust_query)) {
    gene_i <- top_de$gene[top_de$cluster==clust_query[i]]
    ovlp_res <- sapply(marker_list, function(marker_x) {
      ovlp_test <- GeneOverlap::newGeneOverlap(listA = gene_i, listB = marker_x, genome.size = length(unique(bk_gene))) %>% 
        GeneOverlap::testGeneOverlap()
      # if (is.infinite(fshexact_test$estimate)) fshexact_test <- fisher.test(matrix(c(x11, x21, x12, x22)+1, nrow = 2, ncol = 2))
      c(pval=ovlp_test@pval, or=ovlp_test@odds.ratio)
    })
    clust_pval_df[, i] <- ovlp_res[1, rownames(clust_pval_df)]
    clust_or_df[, i] <- ovlp_res[2, rownames(clust_pval_df)]
  }
  return(list(pval=clust_pval_df, or=clust_or_df))
}
