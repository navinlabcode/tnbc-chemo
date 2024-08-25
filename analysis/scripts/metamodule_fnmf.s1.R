# Basic QC of NMF factors
# Author: Yun Yan
library(tidyverse)
library(readr)

tmp <- Sys.glob(file.path(
    '.', 'metamodule_fnmf',
    'nmf_each_sample', 'ARTC*', 'K*', 'reconstruction_err.txt'))
print(length(tmp))

for (i in seq_along(tmp)) {
    dir_res <- dirname(tmp[i])
    f_cell2nmf  <- file.path(dir_res, 'deliver.cell_best_nmf.rds')
    f_nmfmarker <- file.path(dir_res, 'deliver.nmf_markers.rds')
    f_nmftopgene <- file.path(dir_res, 'deliver.nmf_markers.top.rds')
    
    cell2nmf <- read_rds(f_cell2nmf)
    df1 <- enframe(c(table(cell2nmf)), name = 'fnmf', value = 'nCell')
    
    nmfmarker <- read_rds(f_nmfmarker); #str(nmfmarker)
    df2 <- enframe(sapply(nmfmarker, length), name='fnmf', value='nMarkerGene')
    df <- left_join(df1, df2, by='fnmf')
    
    df$nCell[is.na(df$nCell)] <- 0
    df$nMarkerGene[is.na(df$nMarkerGene)] <- 0
    use <- df$nCell >= 2 & df$nMarkerGene >= 3 # beta-version [in use]
    sum(use)
    res <- df$fnmf[use]
    
    write_lines(res, file.path(dir_res, 'nmf_names_in_use.txt'))
    write_csv(df, file.path(dir_res, 'nmf_metainfo.csv') )
}
