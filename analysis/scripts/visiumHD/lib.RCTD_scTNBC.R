suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal:Create a library of RCTD for TNBC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025-05-13
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
suppressPackageStartupMessages({
    library(readr)
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    theme_set(theme_pubr(base_size = 6, legend = "right") %+replace% theme(axis.ticks.length = unit(0.1, "inch")))
    library(cli)
    library(tictoc)
    library(glue)
    library(scales)
    library(tools)
    library(fs)
    library(nanoparquet)
    library(Seurat)
    library(spacexr)
    library(Rfast)
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {

} else {

}

f_sr3 <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/deliver/atlas_downto1000cells_percellstate/ready.sr3.rds"
sr3 <- read_rds(f_sr3)
print(sr3) # 54,915 cells
print(table(sr3$celltypes))
# Tumor   Mye     T     B Fibro  Endo  Peri  <NA> 
# 11357 13676 14000  5330  4000  3872  2680     0 
set.seed(1024)
cells <- intersect(sample(Cells(sr3),size = 10000, replace = F), Cells(sr3))
sr3_10k <- subset(sr3, cells = cells)
print(sr3_10k)
print(table(sr3_10k$celltypes))
#     B  Endo Fibro   Mye  Peri     T Tumor 
#   943   668   716  2523   516  2555  2079 
## build reference
counts <- GetAssayData(sr3_10k, assay = "RNA", slot = "counts")
cluster <- as.character(sr3_10k$celltypes)
cluster <- factor(cluster, levels=c('Tumor', 'Mye', 'T', 'B', 'Fibro', 'Endo', 'Peri'))
names(cluster)  <- Cells(sr3_10k)
print(head(cluster))
print(table(cluster, useNA = "always"))

nUMI <- sr3_10k$nCount_RNA
head(nUMI)
print(range(nUMI))
# 500 222367

length(nUMI) == length(cluster)
ncol(counts) == length(cluster)

# RCTD only max accepts number of cells to: 10000
reference <- Reference(counts, cluster, nUMI)

saveRDS(
    reference, 
    file = "/volumes/USR1/yyan/project/tnbc_visium_hd/lib/RCTD_Lib/RCTD_object.TNBC_celltypes.rds")




cat("[done]")
timestamp()
