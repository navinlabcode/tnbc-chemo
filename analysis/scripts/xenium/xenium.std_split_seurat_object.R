#---------------------------
# Split seurat objects              ----    
#---------------------------
library(Seurat)
library(tidyverse)
library(readr)
options <- commandArgs(trailingOnly = TRUE)
source('/volumes/USR/yyan7/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R')

# f_in <- '/volumes/USR/yyan7/project/tnbc_xenium/data/ART23/nonbinarized_pca/xenium_nonbinarized_pca.seurat.rds'
# f_cluster <- '/volumes/USR/yyan7/project/tnbc_xenium/data/ART23/nonbinarized_pca/celltyping.transferlabel/df.TransferData.Diagnosis.celltype.rds'
# by_what <- 'predicted.id.diagnosed'
f_in <- by_what <- f_cluster <- NULL
if (length(options)==2) {
  f_in  <- options[[1]]
  by_what <- options[[2]]
  dir_res <- file.path(
    dirname(f_in),
    sprintf('objects_split_into_%s', by_what))
}

if (length(options)==3) {
  f_in  <- options[[1]]
  by_what <- options[[2]]
  f_cluster <- options[[3]]
  dir_res <- file.path(
    dirname(f_cluster),
    sprintf('objects_split_into_%s', by_what))
}


fs::dir_create(dir_res)

sr3_giant <- read_rds(f_in); print(sr3_giant)

if (!is.null(f_cluster)) {
  df_cluster <- read_rds(f_cluster)
  stopifnot(all(rownames(df_cluster) %in% Cells(sr3_giant)))
  identical( rownames(df_cluster), Cells(sr3_giant) )
  sr3_giant <- AddMetaData(sr3_giant, df_cluster)
}


Idents(sr3_giant) <- by_what

print(table(Idents(sr3_giant), useNA = 'ifany'))
print(levels(Idents(sr3_giant)))

instance_vector <- levels(Idents(sr3_giant))
ncells_xx <- rep(0, length(instance_vector)); names(ncells_xx) <- instance_vector

for (xx in instance_vector) {
  cat(xx, '')
  dir_saveto <- file.path(dir_res, xx)
  saveto <- file.path(dir_saveto, 'pass.xenium.seurat.rds')
  fs::dir_create(dir_saveto)
  
  idx <- c(sr3_giant[[]][, by_what] %in% xx)
  ncells_xx[xx] <- sum(idx, na.rm = T)
  
  if (file.exists(saveto)){
    cat('existed...')
    next()
  }
  
  sr3 <- subset(sr3_giant, cells = Cells(sr3_giant)[idx])
  cat(ncol(sr3), 'cells...')
  sr3@project.name <- xx
  
  fs::dir_create(dir_saveto)
  cat('export...')
  write_xenium_rds(sr3, dir_saveto, obj_type = 'pass')
  
  rm(sr3)
  
  cat('[done]\n')
}


if ( all(ncells_xx != table(sr3_giant[[by_what]])) ) {
  warning('Cell num does NOT match. Recheck split!')
}

enframe(ncells_xx, name='ident', value='ncells') %>%
  write_csv(file.path(dir_res, 'ncells.csv'))
