# Run CopyKat on ST
#
# Yiyun Lin
#
#


#this is not run in singularity, run in 8787 R4.1.2
rm(list = ls())
library(copykit)
library(readr)
library(base)
library(reader)
library(stringr)
library(dplyr)
library(magrittr)
library(BiocParallel)
source("/volumes/USR1/yiyun/Script/IndexART.R")
source("/volumes/USR1/yiyun/Script/scRNA/scRNA_packages.R")

seu_objlist = read_rds('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_seuratRDS/object_list.rds') 

cpk_objlist = list()

ref = '_nonromal' 
all_samples = c('ART65','ART314','ART31','ART43','ART18','ART40','ART10','ART23','ART304','ART305','ART311','ART312')
sample_inpaper = c('ART65','ART10','ART23','ART304','ART305','ART311','ART312')
all_samples[all_samples %!in% sample_inpaper]
for(pt.name in all_samples){
  samplename=paste0(pt.name,ref)
  print(samplename)
  #read data
  if(T){    
    wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/copycat/',samplename,'/'); dir.create(wkdir); setwd(wkdir)
    tumor <- readVarbinCNA_copykat(paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/copycat/',samplename,'/'),remove_Y = T,genome_version = 'hg19',method = "none")
  }
  
  col = list(ref = c('HBCA' = '#C0D6DF', 'Sample' = '#EE6C4D'),
             Tumor_score=c('Tumor_high'='#F79256','Tumor_low'='#7DCFB6','HBCA' = '#C0D6DF'),
             Copykat.pred = c('diploid'= '#009DAE','nondiploid' = '#F58840'),
             inferTumor = c('Tumor'='#F58840','Normal' = '#009DAE'))
  
  #add annotation
  if(T){
    meta = SummarizedExperiment::colData(tumor) %>% data.frame()
    meta$ref = ifelse(grepl('hbca_',meta$sample),'HBCA','Sample')
    meta$ref %>%table()
    dim(meta)
    SummarizedExperiment::colData(tumor)$ref = meta$ref
    
    SeuratObject = seu_objlist[[pt.name]]
    cell_score = SeuratObject@assays$predictions@data%>%t() %>% data.frame()
    cell_score$tumor %>% hist()
    cell_score$cellname = rownames(cell_score)
    
    Tumor_high = cell_score %>% dplyr::filter(tumor>0.1) %>% pull(cellname) %>% tolower()
    Tumor_high = gsub('-','_',Tumor_high)
    meta$Tumor_score= NA
    meta$Tumor_score[meta$sample%in%Tumor_high] = 'Tumor_high'
    meta$Tumor_score[meta$ref == 'HBCA'] = 'HBCA'
    meta$Tumor_score[is.na(meta$Tumor_score)] = 'Tumor_low'
    SummarizedExperiment::colData(tumor)$Tumor_score = meta$Tumor_score 
    
    
    ART304_pred = read.table(paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/copycat/',pt.name,'PreTxCore_prediction',ref,'.txt'),header = T)
    nonD = ART304_pred%>% filter(CopyKat.pred =='nondiploid') %>% pull(cell.names) %>% tolower()
    nonD = gsub('-','_',nonD)
    
    meta$Copykat.pred = NA
    meta$Copykat.pred[meta$sample%in%nonD]= "nondiploid"
    meta$Copykat.pred[meta$sample%!in%nonD]= "diploid"
    SummarizedExperiment::colData(tumor)$Copykat.pred = meta$Copykat.pred 
    meta$Copykat.pred %>% table()
    
    
    #add annotation
    svdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/plots/',pt.name,'/');dir.create(svdir);
    label =c('ref','Tumor_score','Copykat.pred')
    p1 = plotHeatmap_copykat(tumor, label = label,label_colors = col[label],order_cells = 'hclust',n_threads = 100)
    pdf(paste0(svdir,'heatmap_with_normal.pdf'),width = 10,height = 10)
    print(p1)
    dev.off()
  }
  
  #run umap
  if(T){
    tumor <- runUmap(tumor, n_neighbors = 30, min_dist = 0,seed = 17)
    plotUmap(tumor)
    tumor <- findSuggestedK2(tumor,method = 'leiden',B = 200)
    wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS'); dir.create(wkdir); setwd(wkdir)
        
    plotSuggestedK(tumor)
    k = tumor@metadata$suggestedK
    k = 30 #n neighber = 10
    if(samplename %in% c('ART31_nonromal','ART18_nonromal')){
      k =5
    }
    tumor <- findClusters(tumor,
                          method = 'leiden',
                          k_superclones = 30,
                          k_subclones = k)
    
    col = list(ref = c('HBCA' = '#C0D6DF', 'Sample' = '#EE6C4D'),
               subclones = colormaker(colData(tumor)$subclones),
               Tumor_score=c('Tumor_high'='#F79256','Tumor_low'='#7DCFB6','HBCA' = '#C0D6DF'),
               Copykat.pred = c('diploid'= '#009DAE','nondiploid' = '#F58840'),
               inferTumor = c('Tumor'='#F58840','Normal' = '#009DAE'))
    
    # tumor <- tumor[,colData(tumor)$subclones != 'c0']
    lvls = levels(colData(tumor)$subclones)[levels(colData(tumor)$subclones) %!in% 'c0']
    
    colData(tumor)$subclones = factor(colData(tumor)$subclones, levels = lvls)
    plotUmap(tumor,label = c('subclones'))+
      scale_fill_manual(values = col$subclones,breaks = names(col$subclones))+  
      labs(fill = "subclones") 
    
    plotUmap(tumor,label = c('ref'))+
      scale_fill_manual(values = col$ref,breaks = names(col$ref))+  
      labs(fill = "ref") 
    
    write_rds(tumor,paste0(samplename,'_cpk.rds'))
    
    label = c('ref','subclones','Tumor_score','Copykat.pred')
    p2 = plotHeatmap_copykat(tumor, label = label, label_colors = col[label],row_split = 'subclones',n_threads = 100)
    pdf(paste0(svdir,'heatmap_with_normal_subclone.pdf'),width = 10,height = 10)
    print(p2)
    dev.off()
  }
  
  #distinguish if it's a normal
  if(T){
    wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS'); dir.create(wkdir); setwd(wkdir)
    tumor = read_rds(paste0(samplename,'_cpk.rds'))
    
    tumor <- calcConsensus(tumor,consensus_by = 'subclones')
    cluster_consensus = copykit::consensus(tumor) %>% t() %>% data.frame();cluster_consensus[,1:5]
    cluster_consensus[1:5,1:5]
    seg_data <- t(SummarizedExperiment::assay(tumor, 'segment_ratios'))
    seg_data[1:5,1:5]  
    seg_data = cluster_consensus
    consensus_mean = cluster_consensus %>% colMeans()
    #calculate the mean of square
    msq = apply(seg_data,MARGIN = 1,FUN = function(x){a = sum(x^2)/ncol(seg_data);return(a)})
    msq %>% hist()
    #calculate the correlation
    corr = apply(seg_data,MARGIN = 1, FUN = function(x){
      a = cor.test(x,unname(consensus_mean))
      b = a$estimate %>% unname()
      return(b)})
    corr %>% hist()
    df = cbind(msq,corr) %>% data.frame()%>% mutate(subclones = rownames(cluster_consensus))
    
    
    non_di1 = df$msq > 0.004
    if(pt.name %in% c('ART312','ART40')){
      non_di1 = df$msq > 0.004
    }
    if(ref == '_nonromal'){
      non_di2 = df$corr > 0.3
    }else{
      non_di2 = df$corr > 0.5
    }
    df$inferTumor = NA
    df$inferTumor[non_di1&non_di2] = 'Tumor'
    df$inferTumor[is.na(df$inferTumor)] = 'Normal'
    df$inferTumor = factor(df$inferTumor, levels = c('Tumor','Normal'))

    #check which clone has ref normal
    meta = SummarizedExperiment::colData(tumor) %>% as.data.frame()
    meta2 = left_join(meta,df)
    SummarizedExperiment::colData(tumor)$msq_scl = meta2$msq
    SummarizedExperiment::colData(tumor)$corr_scl = meta2$corr
    SummarizedExperiment::colData(tumor)$inferTumor = meta2$inferTumor
    table(meta2$inferTumor) %>% print()
    
    label = c('ref','Tumor_score','Copykat.pred','inferTumor')
    p3 = plotHeatmap_copykat(tumor, label = label, label_colors = col[label], row_split = 'inferTumor',n_threads = 100)
    pdf(paste0(svdir,'heatmap_with_normal_inferTumor.pdf'),width = 10,height = 10)
    print(p3)
    dev.off()
  }
  
  #save the result in object list
  cpk_objlist[[pt.name]] = tumor
  
  wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/'); dir.create(wkdir); setwd(wkdir)
  write_rds(cpk_objlist,'object_list.rds')
}

###check! if tumor definition same as before
if(T){
  wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/'); dir.create(wkdir); setwd(wkdir)
  cpk_objlist = read_rds('object_list.rds')
  for(pt.name in all_samples){
    samplename=paste0(pt.name,ref)
    message(samplename)
    meta_all = read_rds('/volumes/USR1/yiyun/Project/ARTEMIS_ST/Figure_for_paper/copykat/all_tumor_metadata.rds')
    
    message('old')
    meta_all %>% filter(patient == pt.name) %>% pull(inferTumor) %>% table() %>% print()
   
    message('new')
    tumor = cpk_objlist[[pt.name]]
    meta = colData(tumor) %>% data.frame()
    meta$inferTumor %>% table() %>% print()
    
    message('overlap')
    a = meta_all %>% filter(patient == pt.name,inferTumor == 'Tumor') %>% pull(sample)
    b = meta %>% filter(inferTumor == 'Tumor') %>% pull(sample)
    sum(a %in% b) %>% print()
  }
}




