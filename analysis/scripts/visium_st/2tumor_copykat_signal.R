## Calcualte MSQ score for the inferred CNA of Visium
##
## Yiyun Lin
##
if(T){
  options(stringsAsFactors = F)
  source("scRNA_packages.R")
  library("SChartPack")
  library(copykit)
  library("akima")
  library(patchwork)
  library("randomForestSRC")
  library("packcircles")
  library("dplyr")
  library(readr)
  library("magrittr")
  library("dbscan")
  library("pheatmap")
  library("spatstat")
  library("Seurat",lib.loc="/volumes/USR1/yiyun/R/x86_64-redhat-linux-gnu-library/4.0") 
  # library("SeuratData")
  library("reshape2")
  library("visNetwork")
  library("shiny")
  library("plotly")
  library("viridis")
  library("RColorBrewer")
  # library("ConsensusClusterPlus")
  library("philentropy")
  library('readr')
  `%!in%` = Negate(`%in%`)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(tidyr)
  library(rstatix)
  library(ggpubr)
  library(ggprism)
  # library('hdf5r',lib.loc = '/volumes/USR3/minhu/anaconda3/envs/seurat4/lib/R/library/')
}




#calculate the msq, corr for every patient by spot, use msq>=0.045 as cutoff
wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/'); dir.create(wkdir); setwd(wkdir)
cpk_objlist = read_rds('object_list.rds')
if(TRUE){
  for(pt.name in c('ART65','ART314','ART31','ART43','ART18','ART40','ART10','ART23','ART304','ART305','ART311','ART312')){
    reff ='_nonromal'
    
    samplename=paste0(pt.name,reff)
    svdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/plots/',pt.name,'/');dir.create(svdir);
    print(samplename)
    
    #read
    tumor = cpk_objlist[[pt.name]]
    if(T){
      colData(tumor)$msq = NULL
      colData(tumor)$corr = NULL
      colData(tumor)$patient = NULL
      colData(tumor)$confidence = NULL
    }
    
    #reannotate ART31 and ART18 since they has no normal
    if(pt.name %in% c('ART31', 'ART18')){
      colData(tumor)$inferTumor = 'Normal'
    }
    meta = colData(tumor) %>% data.frame();head(meta)
    if(sum(colData(tumor)$inferTumor == 'Tumor')>0){
      #gate out tumor for further clustering because we want high purity tumor spots
      #for tumor
      if(T){
        tumor_sub <- tumor[,colData(tumor)$inferTumor == 'Tumor']
        tumor_sub <- runUmap(tumor_sub, n_neighbors = 30, min_dist = 0,seed = 17)
        k=30
        tumor_sub <- findClusters(tumor_sub,
                                  method = 'leiden',
                                  k_superclones = 30,
                                  k_subclones = k)
        colData(tumor_sub)$subclones %>% table()%>% print()
      }
      
      #for normal
      if(T){
        normal_sub <- tumor[,colData(tumor)$inferTumor == 'Normal']
        colData(normal_sub)$subclones = 'Normal'
      }
      
      #merge tumor normal again, and calculate msq, corr for each spot
      if(T){
        tumor_up = cbind(tumor_sub,normal_sub)
        seg_data <- t(SummarizedExperiment::assay(tumor_up, 'segment_ratios'))
        seg_data[1:5,1:5]  
        consensus_mean = seg_data %>% colMeans()

        #calculate the mean of square for each spot
        msq = apply(seg_data,MARGIN = 1,FUN = function(x){a = sum(x^2)/ncol(seg_data);return(a)})
        msq %>% hist()
        
        #calculate the correlation for each spot
        corr = apply(seg_data,MARGIN = 1, FUN = function(x){
          a = cor.test(x,unname(consensus_mean))
          b = a$estimate %>% unname()
          return(b)})
        corr %>% hist()
        
        df = cbind(msq,corr) %>% data.frame()
        df$sample = names(msq)
        
        meta = colData(tumor_up) %>% data.frame()
        meta = left_join(meta,df) %>% mutate(patient = pt.name)
        df_cut = meta %>% dplyr::group_by(subclones) %>% dplyr::summarise(msq = mean(msq))
        df_cut$confidence = NA
        df_cut$confidence[df_cut$msq >= 0.045] ='High_purity'
        df_cut$confidence[df_cut$msq < 0.045 & df_cut$msq >= 0.020] ='Medium_purity'
        df_cut$confidence[df_cut$msq < 0.020] ='Low_purity'
        df_cut$confidence[df_cut$subclones == 'Normal'] = 'Normal'
        meta = left_join(meta,df_cut %>% select(-msq))
        meta$confidence = factor(meta$confidence, levels = c('Normal','Low_purity','Medium_purity','High_purity'))
        meta$confidence %>% table() %>% print()
        
        
        colData(tumor_up) = cbind(colData(tumor_up), meta[,colnames(meta) %!in% colnames(colData(tumor_up))])
        label = c('subclones','Tumor_score','Copykat.pred','inferTumor','confidence')
        
        
        col = list(ref = c('HBCA' = '#C0D6DF', 'Sample' = '#EE6C4D'),
                   subclones = colormaker(colData(tumor_up)$subclones),
                   Tumor_score=c('Tumor_high'='#F79256','Tumor_low'='#7DCFB6','HBCA' = '#C0D6DF'),
                   Copykat.pred = c('diploid'= '#009DAE','nondiploid' = '#F58840'),
                   inferTumor = c('Tumor'='#F58840','Normal' = '#009DAE'),
                   area = c('a1' = '#D9534F','a2' = '#FFAD60', 'a3' = '#FFEEAD','a4'='#FFE162','a5' = '#FFFDDE','Normal' = '#009DAE'),
                   confidence = c('High_purity' = '#D9534F','Medium_purity' = '#FFAD60','Low_purity' = '#FFEEAD','Normal' = '#009DAE'))
        
        p1 = plotHeatmap_copykat(tumor_up, label = label, label_colors = col[label],row_split = 'subclones',n_threads = 100)
        pdf(paste0(svdir,'heatmap_subclone_based_on_tumor_only.pdf'),width = 10,height = 10)
        print(p1)
        dev.off() 
        p2 = plotHeatmap_copykat(tumor_up, label = label, label_colors = col[label],row_split = 'confidence',n_threads = 100)
        pdf(paste0(svdir,'heatmap_confidence.pdf'),width = 10,height = 10)
        print(p2)
        dev.off() 
      }
    }else{
      tumor_up = tumor
    }
    cpk_objlist[[pt.name]] = tumor_up
  }
}
wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/'); dir.create(wkdir); setwd(wkdir)
write_rds(cpk_objlist,'object_list.rds')

meta_all = data.frame()
for(pt.name in c('ART10','ART23','ART43','ART65','ART304','ART305','ART311','ART312','ART314','ART40','ART31','ART18')){
  tumor = cpk_objlist[[pt.name]]
  if(pt.name %in% c('ART31','ART18')){
    colData(tumor)$subclones = 'Normal'
    colData(tumor)$msq = 0
    colData(tumor)$corr = 0
    colData(tumor)$patient = pt.name
    colData(tumor)$confidence = 'Normal'
  }
  cpk_objlist[[pt.name]] = tumor
  meta = colData(tumor) %>% data.frame()
  meta_all = rbind(meta_all,meta)
}
wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/'); dir.create(wkdir); setwd(wkdir)
write_rds(cpk_objlist,'object_list.rds')
write_rds(meta_all,'all_tumor_metadata.rds')


meta_all_new = read_rds('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/all_tumor_metadata.rds')
table(meta_all_new$patient,meta_all_new$subclones)

# Justify why msq>0.045 is used to select high-tumor purity spots
# important note, area is just subclone re-ordered by msq score
pdf('/volumes/USR1/yiyun/Project/ARTEMIS_ST/pass_to_Yun_05152023/origin_copykitRDS/msq_box_for_cutoff_selection.pdf',width = 15,height = 5)
ggboxplot(meta_all_new %>% dplyr::filter(subclones != 'Normal'),x = 'subclones',y='msq') + 
  facet_grid(~patient) + 
  geom_hline(yintercept = 0.045)+ 
  geom_hline(yintercept = 0.020)
dev.off()

if(TRUE){
  wkdir = paste0('/volumes/USR1/yiyun/Project/ARTEMIS_ST/Figure_for_paper/copykat/'); dir.create(wkdir); setwd(wkdir)
  meta_all = read_rds(file = 'all_tumor_metadata.rds')
  meta_all$area = factor(meta_all$area,levels = c(paste0('a',1:5),'Normal'))
  #pick cutoff according to quantile
  meta_all %>% filter(inferTumor  == 'Tumor') %>% pull(msq) %>% quantile()
    
  df_cut = meta_all %>% dplyr::group_by(patient,area) %>% dplyr::summarise(msq = median(msq))
  
  df_cut$confidence = df_cut$area %>% varhandle::unfactor()
  df_cut$confidence[df_cut$msq >= 0.045] ='High_purity'
  df_cut$confidence[df_cut$msq < 0.045 & df_cut$msq >= 0.020] ='Medium_purity'
  df_cut$confidence[df_cut$msq < 0.020] ='Low_purity'
  df_cut$confidence[df_cut$area == 'Normal'] = 'Normal'
  table(df_cut$confidence)
  
  meta_all  = left_join(meta_all, df_cut %>% dplyr::select(-msq))
  meta_all$confidence = factor(meta_all$confidence,levels = c('High_purity','Medium_purity','Low_purity','Normal'))
  table(meta_all$confidence)
  
  # meta_all %>% write_rds(., file = 'all_tumor_metadata.rds')
}
