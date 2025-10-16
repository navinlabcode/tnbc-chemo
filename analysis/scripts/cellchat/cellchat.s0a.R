#---------------------------
# Downsampling the original dataset              ----  
# - Avoid patient bias
# - Enough cell number per cell state
#---------------------------
library(Seurat)
library(tidyverse)

#------------------- ~~~ Downsampling ~~~ -------------------  

f_in <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/ready.sr3.rds'
dir_res <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/ecotype/cellchat'
fs::dir_create(dir_res)
dict_update_cellstatenames <- read_tsv('/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/color.css', col_names = c('old', 'color', 'new'), col_types = c('c', 'c', 'c'))
dict_update_cellstatenames <- dict_update_cellstatenames %>% filter(!new %in% c('Unresolved', NA) )
cellstatenames_old <- dict_update_cellstatenames$old # 60 features
cellstatenames_new <- dict_update_cellstatenames$new # 60 features

sr3 <- read_rds(f_in)
print(sr3) # 427857 cells

table(sr3$ecotrait_feature)
sum(table(sr3$patient, useNA='ifany'))


mat_table <- table(sr3$patient, sr3$ecotrait_feature) %>% as.data.frame.matrix()
mat_table <- mat_table[, match(cellstatenames_old, colnames(mat_table))]
mat_table <- as.matrix(mat_table)
view(mat_table)

thre_ncell_per_patient_per_cstate <- ceiling(colMeans(mat_table)) #final yield: 208,508/427,857
thre_ncell_per_patient_per_cstate <- ceiling(colMeans(mat_table)/2) # 132,318/427,857
thre_ncell_per_patient_per_cstate <- ceiling(colMeans(mat_table)/3) # 97,777/427,857
thre_ncell_per_patient_per_cstate <- ceiling(colMeans(mat_table)/5) # 70,011/427,857
thre_ncell_per_patient_per_cstate[thre_ncell_per_patient_per_cstate<10] <- 10

df_cellmeta <- sr3@meta.data

cnames_use <- c()
for (cs in cellstatenames_old) {
  message(cs)
  df_cellmeta_cs <- df_cellmeta[df_cellmeta$ecotrait_feature %in% cs, ]
  cs_cnames <- subset_downsample(
    rownames(df_cellmeta_cs), 
    df_cellmeta_cs$patient, 
    n_at_least = thre_ncell_per_patient_per_cstate[cs]
  )
  message(length(cs_cnames), '/', nrow(df_cellmeta_cs))
  cnames_use <- c(cnames_use, cs_cnames)
}; rm(df_cellmeta_cs)
message(length(cnames_use), '/', nrow(df_cellmeta))

subset_downsample <- function(x, grouping, n_at_least=1) {
  ## Similar to Seurat's subset(..., downsample)
  ## x: vector
  ## grouping: vector
  ## n_at_least: int
  idx <- seq_along(x)
  idx_use <- c()
  
  ident_opts <- unique(grouping)
  for (ident_i in ident_opts) {
    idx_source <- idx[grouping %in% ident_i]
    idx_sample <- NULL
    if (length(idx_source) > n_at_least) {
      idx_sample <- sample(idx_source, size = n_at_least, replace = F)
    } else {
      idx_sample <- idx_source
    }
    idx_use <- c(idx_use, idx_sample)
  }
  return(x[idx_use])
}

df_cellmeta_use <- df_cellmeta[rownames(df_cellmeta) %in% cnames_use, ]
table(df_cellmeta_use$patient, df_cellmeta_use$ecotrait_feature)
sort(table(df_cellmeta_use$patient), decreasing = T)
# ARTC41  ARTC94  ARTC93  ARTC81  ARTC92  ARTC47 
# 966     963     949     938     910     901 
# ARTC77  ARTC62  ARTC72  ARTC73  ARTC97  ARTC40 
# 901     884     876     876     873     867 
# ARTC39  ARTC80  ARTC18  ARTC12  ARTC28  ARTC98 
# 859     857     853     846     843     843 
# ARTC43 ARTC109  ARTC32 ARTC102  ARTC71  ARTC44 
# 837     833     832     829     828     827 
# ARTC84  ARTC27  ARTC35  ARTC87  ARTC95  ARTC67 
# 818     803     803     801     797     794 
# ARTC99  ARTC17  ARTC68  ARTC63 ARTC105  ARTC74 
# 794     791     775     774     773     767 
# ARTC69  ARTC66  ARTC42  ARTC16  ARTC19  ARTC49 
# 764     758     754     751     743     740 
# ARTC24  ARTC10 ARTC101  ARTC65 ARTC106  ARTC89 
# 738     735     734     732     731     731 
# ARTC78  ARTC09  ARTC83  ARTC52  ARTC70 ARTC104 
# 726     724     723     717     716     711 
# ARTC58  ARTC23  ARTC50  ARTC11  ARTC54  ARTC13 
# 705     703     701     690     690     688 
# ARTC85  ARTC61  ARTC91  ARTC79  ARTC88  ARTC90 
# 685     681     678     661     659     655 
# ARTC100  ARTC96  ARTC86 ARTC103 ARTC107  ARTC01 
# 652     640     634     631     629     619 
# ARTC14  ARTC03  ARTC45  ARTC31  ARTC20  ARTC55 
# 619     615     600     598     574     574 
# ARTC48  ARTC75  ARTC08  ARTC82  ARTC26  ARTC29 
# 572     567     543     543     522     521 
# ARTC06  ARTC33  ARTC21  ARTC37  ARTC53  ARTC36 
# 517     515     497     494     494     479 
# ARTC30  ARTC07  ARTC34  ARTC51  ARTC15  ARTC25 
# 469     462     461     442     423     415 
# ARTC38  ARTC05  ARTC76  ARTC04  ARTC02 
# 410     358     353     334     330 
sort(table(df_cellmeta_use$ecotrait_feature), decreasing = T)
# CD8_gzmk            CD4_treg 
# 8597                5092 
# CD4_Th_cxcr4      CD4_naive_ccr7 
# 5016                4008 
# CD8_exh_cxcl13        B_plasma_IgG 
# 3890                3368 
# CD4_exh_cxcl13          NK_cd16low 
# 3154                3137 
# CD4_ifn             B_naive 
# 1756                1525 
# CD8_nklike_xcl1             CD8_ifn 
# 1491                1490 
# B_mem    macro-CXCL10-IFN 
# 1255                1107 
# CD8_Tm_il7r           nkt_temra 
# 1042                 991 
# mye-prol            macro-m1 
# 955                 936 
# pDC         macro-CCL18 
# 929                 921 
# t_prol         NK_cd16high 
# 887                 872 
# M04__Stress          mono-class 
# 862                 837 
# M07__S.G1                cDC2 
# 824                 819 
# M08__Hypoxia     M12__Cholestero 
# 777                 770 
# M13__StressER M10__EpithelialDiff 
# 760                 758 
# M11__LumSec          M09__Basal 
# 755                 753 
# M05__Interferon            M06__HLA 
# 749                 746 
# M01__G2.M        B_plasma_IgA 
# 729                 677 
# cDC1                 mDC 
# 656                 635 
# B_mem_ifn          macro-lipo 
# 571                 539 
# Mast            macro-m2 
# 373                 371 
# endo-Tip         endo-venous 
# 334                 316 
# endo-arterial            peri-myo 
# 298                 283 
# CAFs      endo-capillary 
# 272                 253 
# TAM         peri-immune 
# 251                 217 
# neutro-FCGR3B         peri-matrix 
# 214                 206 
# Fibro       Fibro-PCOLCE2 
# 204                 199 
# Fibro-immuno              B_germ 
# 178                 136 
# endo-prol           peri-crem 
# 109                  81 
# endo-lymph          macro-cd52 
# 40                  40 
sort(table(df_cellmeta$ecotrait_feature), decreasing = T)
table(duplicated(cnames_use))
cnames_use <- rownames(df_cellmeta_use)
str(cnames_use)
dir_res
write_lines(cnames_use, file.path(dir_res, 'cnames_for_LR_analysis.txt'))
write_rds(df_cellmeta_use, file.path(dir_res, 'sr3_metadata.df.rds'))
sr3 <- subset(sr3, cells = cnames_use)
print(sr3)
write_rds(sr3, file.path(dir_res, 'ready.sr3.rds'))

#------------------- ~~~ Export as cellchat object ~~~ -------------------  
data.input <- GetAssayData(sr3, assay = 'RNA', slot = 'data')
v_celltypes <- str_split_i(Cells(sr3), '_', 1)
v_cellstates_old <- as.character(sr3$ecotrait_feature)
v_cellstates_new <- deframe(dict_update_cellstatenames[, c('old', 'new')])[v_cellstates_old]
meta <- data.frame(
  celltype = v_celltypes, 
  cellstate_old = v_cellstates_old,
  cellstate = v_cellstates_new,
  archetype = sr3$archetype, 
  pCR_status = sr3$pCR_status, 
  row.names = Cells(sr3), stringsAsFactors = F)
meta$cellstate <- factor(meta$cellstate, levels = cellstatenames_new)
meta$cellstate_old <- factor(meta$cellstate_old, levels = cellstatenames_old)
meta$group <- meta$cellstate

write_rds(data.input, file.path(dir_res, 'data.mat.rds'))
write_rds(meta, file.path(dir_res, 'sc_metadata.df.rds'))
library(CellChat)

data.input <- read_rds(file.path(dir_res, 'data.mat.rds'))
meta <- read_rds(file.path(dir_res, 'sc_metadata.df.rds'))

#------ load ecotrait and reorder cellstates ------
load(file.path('/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102',
               'ecotype', 'use_tumor_hybrid_frac', 'allpatient_with_archetypes', 'use.RData'))
head(memb_feature)
## manually change and match to the main figure
dict_ecotrait_old2new <- structure(
  c(8, 1, 3, 2, 7, 4, 6, 5),
  names = 1:8)
memb_feature <- structure( 
  ruok::replace_vector(memb_feature, dict_ecotrait_old2new),
  names = names(memb_feature))
memb_feature <- structure(
  paste0('ecotrait', memb_feature),
  names = names(memb_feature))
## update the feature names
dict_cellstate_old2new <- deframe(unique(meta[, c('cellstate_old', 'cellstate')]))
names(memb_feature) <- as.character(dict_cellstate_old2new[names(memb_feature)])
memb_feature <- enframe(memb_feature, 'cellstate', 'ecotrait')
memb_feature <- memb_feature %>% dplyr::arrange(ecotrait, cellstate)
memb_feature <- deframe(memb_feature)
meta$cellstate <- fct_relevel(
  meta$cellstate, 
  levels(reorder(levels(meta$cellstate), memb_feature[levels(meta$cellstate)], unique))
)
cellchat <- createCellChat(object = data.input, meta = meta, group.by='cellstate')
table(cellchat@idents)
write_rds(cellchat, file.path(dir_res, 'cellchat.rds'))



