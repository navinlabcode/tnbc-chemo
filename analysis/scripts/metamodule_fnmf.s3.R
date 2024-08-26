##
## Determine marker genes for metaprograms
##
## Yun Yan

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(gtools)
})

dir_proj <- './metamodule_fnmf/'


if (T) {   
    f_mm_nmf <- file.path(dir_proj, 'all_ranks', 'jaccard', 
                          'manual_alt_clean_wardD2', 
                          'jaccard-manual.metamodule_membership.csv')
    mm_nmf <- read_csv(f_mm_nmf)
    dir_res <- file.path(dir_proj, 'MM_alt_clean_byscore_wardD2')
    dir_create(dir_res)
}


mm_names <- mixedsort( unique(mm_nmf$meta_module) ); mm_names

#------ Create MM ~ gene_loadings ------
if (!file.exists(file.path(dir_res, 'fnmf_w.factedBy_metamodule.list.rds'))) {
    mm_loadings <- lapply(mm_names, function(M){
        df_M <- mm_nmf %>% dplyr::filter(meta_module == M) 
        M_nmf <- df_M$index_name
        M_loadings <- sapply(seq_along(M_nmf), function(i)  {
            s <- df_M$sample[i]
            r <- paste0('K', df_M$Rank[i])
            f <- df_M$NMF[i]
            stopifnot(length(s)==length(r))
            w <- read_rds( file.path(dir_proj, 'nmf_each_sample', 
                                     s, r, 'W.matrix.rds') )
            w[, f, drop=T]
            
        })
        colnames(M_loadings) <- M_nmf
        M_loadings
    })
    names(mm_loadings) <- mm_names
    
    write_rds(mm_loadings, file.path(dir_res, 'fnmf_w.factedBy_metamodule.list.rds'))
} else {
    mm_loadings <- read_rds( file.path(dir_res, 'fnmf_w.factedBy_metamodule.list.rds') )
}

if (! file.exists(file.path(dir_res, 'W.matrix.rds'))) {
    mm_loadings_avg <- sapply(mm_loadings, function(x){
        rowMeans(x)
    })
    write_rds(mm_loadings_avg, file.path(dir_res, 'W.matrix.rds'))
} else {
    mm_loadings_avg <- read_rds(file.path(dir_res, 'W.matrix.rds'))
}


#------ Decide markers per mm ------
source('util.nmf.viz.R')
mm_markers <- define_marker_gene_per_nsnmf(mm_loadings_avg, n_tolerance = 0)
str(mm_markers)
mm_top_genes <- top_marker_gene_per_nsnmf(mm_loadings_avg, 50)
str(mm_top_genes)
mm_markers_tol1 <- define_marker_gene_per_nsnmf(mm_loadings_avg, n_tolerance = 1)
str(mm_markers_tol1)
mm_markers_tol2 <- define_marker_gene_per_nsnmf(mm_loadings_avg, n_tolerance = 2)
str(mm_markers_tol2)

write_rds( mm_markers, file.path(dir_res, 'deliver.mm_markers.rds'))
write_gmx2(mm_markers, file.path(dir_res, 'deliver.mm_markers.csv'))

write_rds( mm_markers_tol1, file.path(dir_res, 'deliver.mm_markers_tol1.rds'))
write_gmx2(mm_markers_tol1, file.path(dir_res, 'deliver.mm_markers_tol1.csv'))

write_rds( mm_markers_tol2, file.path(dir_res, 'deliver.mm_markers_tol2.rds'))
write_gmx2(mm_markers_tol2, file.path(dir_res, 'deliver.mm_markers_tol2.csv'))


write_rds( mm_top_genes, file.path(dir_res, 'deliver.mm_top_genes.rds'))
write_gmx2(mm_top_genes, file.path(dir_res, 'deliver.mm_top_genes.csv'))

if (F) {
    ## deprecated -- t.test does not help
    mm_loadings <- do.call('cbind', mm_loadings)
    head(mm_nmf)
    mm_markers_ttest <- define_marker_gene_per_nmf_by_ttest(
        x=mm_loadings, INDEX = deframe(mm_nmf[, c('index_name', 'meta_module')]))
    write_rds( mm_markers_ttest, file.path(dir_res, 'deliver.mm_markers_ttest.rds'))
    write_gmx2(mm_markers_ttest, file.path(dir_res, 'deliver.mm_markers_ttest.csv'))
    str(mm_markers_ttest)
}




    
