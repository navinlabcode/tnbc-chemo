#
# NMF factor similarity matrix 
# - jaccard based on the marker genes
# - pearson correlation based on the gene loadings
#
# Yun Yan (yun.yan@uth.tmc.edu)
#

suppressPackageStartupMessages({
    library(Matrix)
    library(ggpubr); library(ggplot2)
    library(tidyverse)
    library(ComplexHeatmap)
    library(matrixStats)
    library(scales)
    library(WGCNA)
    theme_set(theme_pubr())
    library(future)
    library(future.apply)
})
dir_res <- './metamodule_fnmf/all_ranks'
fs::dir_create(dir_res)

#------ Glob samples with good NMFs of all ranks ------
dir_proj <- './metamodule_fnmf/'    
tmp <- Sys.glob(file.path(
    dir_proj,
    'nmf_each_sample', 'ARTC*', 'K*', 'nmf_names_in_use.txt'))

sample_name_list <- unique(basename(dirname(dirname(tmp))))
str(sample_name_list) #97

# blacklist of patients
# Patients withdrewing constent are removed from all analysis
# Sample without enough cancer cells are removed from metaprogam analysis
df_pat <- read_rds('./sr3_metadata.df.rds')
ncells_pat <- c(table(df_pat$patient))
pat_with_enough_cells <- names(ncells_pat)[ncells_pat >= 20]; str(pat_with_enough_cells)
df_pat <- unique( df_pat[, c('patient', 'PCR_status')] )
pat_wanted <- df_pat$patient[df_pat$PCR_status!='Removed']
table(df_pat$PCR_status)
str(pat_wanted) # 97

table(sample_name_list %in% unique(intersect(pat_wanted, pat_with_enough_cells)))
sample_name_list <- intersect(sample_name_list, unique(intersect(pat_wanted, pat_with_enough_cells)))
str(sort(sample_name_list)) ## 88 patients to use

#------ Create a giant jaccard similarity matrix (factors x factors) ------
## remove some noise factors
source("util.nmf.viz.R")
all_nmfs_markers <- list() 
for ( i in seq_along(tmp) ) {
    
    xx <- tmp[i]
    xx_dir <- dirname(xx)
    xx_sample <- basename(dirname(dirname(xx)))
    xx_nmf_R <- as.numeric(str_remove(basename(dirname(xx)), 'K'))
    if (! xx_sample %in% sample_name_list) { next() }
    if ( is.na(xx_nmf_R) ) { next() }
    cat(xx_sample, xx_nmf_R, '...')
    
    xx_nmf_names <- read_lines(xx)
    # str(xx_nmf_names)
    nmf_markers <- read_rds(file.path(xx_dir, 'deliver.nmf_markers.rds'))
    nmf_markers <- nmf_markers[xx_nmf_names]
    xx_nmf_id <- sprintf('%s_%s_%s', xx_sample, xx_nmf_R, names(nmf_markers))
    names(nmf_markers) <- xx_nmf_id
    
    all_nmfs_markers <- c(all_nmfs_markers, nmf_markers)
}
rm(nmf_markers)
print(length(all_nmfs_markers))


plan(multisession, workers=30)
if (!file.exists(file.path(dir_res, 'jaccard_mat.all_ranks.rds'))) {
    jaccard_mat <- self_pairwise_run(all_nmfs_markers, func = calc_jaccard, 
                                     run_parallel = T)
    write_rds(jaccard_mat, file.path(dir_res, 'jaccard_mat.all_ranks.rds'))
    
    nmfs_max_jaccard <- apply(jaccard_mat, 1, function(x) max(x[x!=1])) ## second closest factor
    pdf(file.path(dir_res, 'jaccard.qc.histogram.nmf_second_nearest_factor.pdf'), width = 6, height = 5, useDingbats = F)
    hist(nmfs_max_jaccard, 100, xlab = 'Jaccard with the second closest nmf factor', ylab='num of query nmf factor')
    dev.off()
    
    print(quantile(nmfs_max_jaccard))
    is_good <- nmfs_max_jaccard >= .25
    print(table(is_good))
    nmfs_names_good <- rownames(jaccard_mat)[is_good]
    write_lines(nmfs_names_good, file.path(dir_res, 'jaccard.nmfs_examined.txt'))
    
    jaccard_mat <- jaccard_mat[nmfs_names_good, nmfs_names_good]
    jaccard_mat <- as.matrix(jaccard_mat)
    print(dim(jaccard_mat))
    write_rds(jaccard_mat, file.path(dir_res, 'jaccard_mat.all_ranks.in_use.rds'))
    
} else {
    jaccard_mat <- read_rds( file.path(dir_res, 'jaccard_mat.all_ranks.rds') )
}
plan(sequential)    

#------ Create nmf meta info data frame (fnmf ~ nCells ~ nMarkerGenes) ------
if (! file.exists(file.path(dir_res, 'nmf_metainfo.csv'))) {
    all_nmfs_meta <- list() 
    for ( i in seq_along(tmp) ) {
        
        xx <- tmp[i]
        xx_dir <- dirname(xx)
        xx_sample <- basename(dirname(dirname(xx)))
        xx_nmf_R <- as.numeric(str_remove(basename(dirname(xx)), 'K'))
        if (! xx_sample %in% sample_name_list) { next() }
        if ( is.na(xx_nmf_R) ) { next() }
        cat(xx_sample, xx_nmf_R, '...')
        
        xx_nmf_names <- read_lines(xx)
        if (length(xx_nmf_names) == 0) { next() }
        xx_df <- read_csv(file.path(xx_dir, 'nmf_metainfo.csv'))
        xx_df <- xx_df[xx_df$fnmf %in% xx_nmf_names, ]
        xx_nmf_id <- sprintf('%s_%s_%s', xx_sample, xx_nmf_R, xx_df$fnmf)
        xx_df$fnmf_id <- xx_nmf_id
        
        xx_H <- read_rds(file.path(xx_dir, 'H.matrix.rds')); dim(xx_H)
        xx_n_total_cells <- ncol(xx_H)
        xx_df$propCell <- xx_df$nCell / xx_n_total_cells
        
        all_nmfs_meta <- c(all_nmfs_meta, list(xx_df))
        
    }
    all_nmfs_meta <- do.call('rbind', all_nmfs_meta)
    write_csv(all_nmfs_meta, file.path(dir_res, 'nmf_metainfo.csv'))
    write_rds(all_nmfs_meta, file.path(dir_res, 'nmf_metainfo.rds'))
    
}
#------ Create a giant gene loading matrix (genes x all factors) ------
if (! file.exists(file.path(dir_res, 'gene_loadings.mat.rds'))) {
    all_nmfs_gloading <- list() 
    for ( i in seq_along(tmp) ) {
        
        xx <- tmp[i]
        xx_dir <- dirname(xx)
        xx_sample <- basename(dirname(dirname(xx)))
        xx_nmf_R <- as.numeric(str_remove(basename(dirname(xx)), 'K'))
        if (! xx_sample %in% sample_name_list) { next() }
        if ( is.na(xx_nmf_R) ) { next() }
        cat(xx_sample, xx_nmf_R, '...')
        
        xx_nmf_names <- read_lines(xx)
        if (length(xx_nmf_names) == 0) { next() }
        xx_W <- read_rds(file.path(xx_dir, 'W.matrix.rds')); print(dim(xx_W))
        xx_W <- xx_W[, xx_nmf_names, drop=FALSE]
        xx_nmf_id <- sprintf('%s_%s_%s', xx_sample, xx_nmf_R, colnames(xx_W))
        colnames(xx_W) <- xx_nmf_id
        
        all_nmfs_gloading <- c(all_nmfs_gloading, list(xx_W))
    }; rm(xx_W); rm(xx_nmf_names)
    stopifnot(all(sapply(2:length(all_nmfs_gloading), function(i) identical(rownames(all_nmfs_gloading[[1]]), rownames(all_nmfs_gloading[[i]]) )))) # make sure all rownames per matrix is the same
    all_nmfs_gloading <- do.call('cbind', all_nmfs_gloading)
    
    write_rds(all_nmfs_gloading, file.path(dir_res, 'gene_loadings.mat.rds'))
} else {    
    all_nmfs_gloading <- read_rds( file.path(dir_res, 'gene_loadings.mat.rds'))
}    
dim(all_nmfs_gloading)

#------ [optional] Create a giant gene-loading correlation matrix (factors x factors) ------
corr_mat <- WGCNA::cor(
    x=all_nmfs_gloading, 
    method='pearson', nThreads=40)
corr_mat <- as.matrix(corr_mat)
write_rds(corr_mat, file.path(dir_res, 'pearson_mat.all_ranks.rds'))

corr_mat <- read_rds( file.path(dir_res, 'pearson_mat.all_ranks.rds') )
dim(corr_mat)
nmfs_max_corr <- apply(corr_mat, 1, function(x) max(x[x!=1])) ## second closest factor
pdf(file.path(dir_res, 'pearson.qc.histogram.nmf_second_nearest_factor.pdf'), width = 6, height = 5, useDingbats = F)
hist(nmfs_max_corr, 100, xlab = 'Pearson corr with the second closest nmf factor', ylab='num of query nmf factor')
dev.off()

print(quantile(nmfs_max_corr))
is_good <- nmfs_max_corr >= 0.7
print(table(is_good))
nmfs_names_good <- rownames(corr_mat)[is_good]
write_lines(nmfs_names_good, file.path(dir_res, 'pearson.nmfs_examined.txt'))

corr_mat <- corr_mat[nmfs_names_good, nmfs_names_good]
corr_mat <- as.matrix(corr_mat)
print(dim(corr_mat))
write_rds(corr_mat, file.path(dir_res, 'pearson_mat.all_ranks.in_use.rds'))
