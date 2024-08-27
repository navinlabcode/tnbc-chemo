## Archetype analysis at pseudo-bulk (patient) level
##
## - identifying archetype
## - designate archetyeps to patients
## - archetype marker genes
## - justify the number of archetypes
## - gene visualization
## - DE analysis
##



library(colorspace)
library(RcppML)
library(Biobase)
library(matrixStats)
library(tictoc); library(cli)
library(readr); library(tidyverse)
library(fgsea)
library(ggpubr); library(ggplot2); library(colorspace); library(ComplexHeatmap)
if (T) {
    ## helper functions
    tau_itai <- function(x) {
        if (any(x<0)) {stop('all values should be non-negative')}
        n <- length(x)
        x_hat <- x / max(x, na.rm=TRUE)
        return(sum(1 - x_hat) / (n-1))
    }
    # tau_itai(c(0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0)) # 0.95
    ratio_of_top2 <- function(x) {
        if (any(x<0)) {stop('all values should be non-negative')}
        a <- max(x, na.rm=TRUE)
        b <- max(x[x!=a], na.rm=TRUE)
        return(a/b)
    }
}
pal_nmf4 <- c(
    'fNMF1' = '#0066CC',
    'fNMF2' = '#99cc00',
    'fNMF3' = '#ff9933',
    'fNMF4' = '#ff00cc'
); show_col(pal_nmf4)

dir_res <- file.path('./psbulk/fastnmf')
fs::dir_create(dir_res)
f_expr <- './psbulk/normalized.matrix.rds'
f_expr_deseq2 <- './psbulk/normalized.DESeq2VST.matrix.rds'
f_in <- f_expr_deseq2
# f_in <- f_expr
mat <- read_rds(f_in)
dim(mat)
#------ Preprare nmf input matrix ------

if (! file.exists(file.path(dir_res, 'input.matrix.rds'))) {
    #------ filter out genes ------
    g_avg <- rowMeans(mat); quantile(g_avg)
    g_var <- rowVars(mat); quantile(g_var)
    g_cv2 <- g_var / (g_avg^2); quantile(g_cv2)
    if ( TRUE ) {
        # idx <- g_avg > quantile(g_avg, .5) & g_var > quantile(g_var, .5) ; print(table(idx))
        idx <- g_avg > quantile(g_avg, .25) & g_var > quantile(g_var, .25) ; print(table(idx))
        idx <- g_avg > quantile(g_avg, .25) & g_avg < quantile(g_avg, .99) & g_var > quantile(g_var, .25) & g_var < quantile(g_var, .99); print(table(idx))
        idx <- g_avg > mean(g_avg) & g_var > quantile(g_var, 0.25) ; print(table(idx))
        idx <- g_avg > mean(g_avg) & g_var > mean(g_var) ; print(table(idx))
    }
    plot(g_avg, sqrt(g_var), col=as.numeric(idx) + 1)
    plot(g_avg, sqrt(g_var))
    points(g_avg[which(names(g_avg) %in% c('CD74'))], 
           sqrt(g_var)[which(names(g_avg) %in% c('CD74'))], col=2, pch=16)
    points(g_avg[which(names(g_avg) %in% c('VCAN'))], 
           sqrt(g_var)[which(names(g_avg) %in% c('VCAN'))], col=3, pch=16)
    mat <- mat[idx, ]
    #------ center and remove negatives ------
    mat <- t(mat)
    mat <- scale(x=mat, center = TRUE, scale = FALSE)
    mat <- t(mat)
    mat[mat < 0] <- 0
    write_rds(mat, file.path(dir_res, 'input.matrix.rds'))
} else {
    mat <- read_rds(file.path(dir_res, 'input.matrix.rds'))
}
print(dim(mat))

#------ Run NMF ------
    
nmf_rank_vec <- c(2:10)
for (r in nmf_rank_vec) {
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    fs::dir_create(tmp)
    if (file.exists(file.path(tmp, 'fastNMF_result_object.rds'))) {
        res <- read_rds(file.path(tmp, 'fastNMF_result_object.rds'))
    } else {
        res <- RcppML::nmf(mat, k=r, tol=1e-5, L1 = c(0.05, 0.05))
        write_rds( res, file.path(tmp, 'fastNMF_result_object.rds') )
    }
    dim(res$w) # genes x factor
    dim(res$h) # factor x samples
    W <- res$w
    H <- res$h
    # W[1:3, 1:3]
    # H[1:3, 1:3]
    
    rownames(W) <- rownames(mat)
    colnames(H) <- colnames(mat)
    colnames(W) <- rownames(H) <- paste0('fNMF', 1:r)
    
    write_rds(W, file.path(tmp, 'W.matrix.rds'))
    write_rds(H, file.path(tmp, 'H.matrix.rds'))
    
}

source('util.ruok.R')
library(ComplexHeatmap)
library(ggpubr)
library(ggplot2)
#------ Decide the dictionary sample ~ NMF ------
for (r in c(4, 5)) {
    cat(r, '...')
    if (r == 5) {
        pal_nmf4 <- c(pal_nmf4, 'fNMF5'='brown')
    }
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    # W <- read_rds(file.path(tmp, 'W.matrix.rds'))
    H <- read_rds(file.path(tmp, 'H.matrix.rds')) # factors x sample
    pat2nmf <- rownames(H)[apply(H, 2, nnet::which.is.max)]; names(pat2nmf) <- colnames(H)
    
    if (F) {
        H_viz <- scale(t(H)) # sample x factor  
        pat2nmf <- colnames(H_viz)[apply(H_viz, 1, which.max)] ## no big differernce from without scale
        names(pat2nmf) <- rownames(H_viz)
    }
    print(table(pat2nmf))
    
    pat2nmf_specificity_tau <- apply(H, 2, tau_itai)
    stopifnot(identical(names(pat2nmf_specificity_tau), names(pat2nmf)))
    pdf(file.path(tmp, 'pat2nmf_specificity.tau.%03d.pdf'), 
        width = 4, height = 4, useDingbats = F, onefile = F)
    vioplot::vioplot(pat2nmf_specificity_tau,ylab='tau')
    boxplot(pat2nmf_specificity_tau ~ pat2nmf, ylab='tau', col = pal_nmf4)
    abline(h=0.6, lty='dashed', col='red')
    beeswarm(pat2nmf_specificity_tau ~ pat2nmf, ylab='tau', col = pal_nmf4, corral = "wrap")
    abline(h=0.6, lty='dashed', col='red')
    dev.off()
    pat2nmf_specificity_fctop2 <- apply(H, 2, ratio_of_top2)
    pat2nmf_specificity_fctop2[is.infinite(pat2nmf_specificity_fctop2)] <- max(pat2nmf_specificity_fctop2[!is.infinite(pat2nmf_specificity_fctop2)])
    pat2nmf_specificity_fctop2 <- log2(pat2nmf_specificity_fctop2)
    stopifnot(identical(names(pat2nmf_specificity_fctop2), names(pat2nmf)))
    pdf(file.path(tmp, 'pat2nmf_specificity.log2_fc_top2.%03d.pdf'), 
        width = 4, height = 4, useDingbats = F, onefile = F)
    vioplot::vioplot(pat2nmf_specificity_fctop2, ylab='log2_fc_top2')
    boxplot(pat2nmf_specificity_fctop2 ~ pat2nmf, ylab='log2_fc_top2')
    dev.off()
    
    df_pat2nmf <- data.frame(fNMF=pat2nmf, 
                             tau=pat2nmf_specificity_tau, 
                             `tau_gt_0.6`=pat2nmf_specificity_tau>0.6)
    row_i <- rownames(df_pat2nmf)[with(df_pat2nmf, order(fNMF, tau))]
    df_pat2nmf <- df_pat2nmf[row_i, ]
    
    row_anno <- HeatmapAnnotation(df=df_pat2nmf, which='row', 
                                  col=list(fNMF=pal_nmf4, 
                                           tau=circlize::colorRamp2(
                                               seq(from=quantile(pat2nmf_specificity_tau, 0.01), to=quantile(pat2nmf_specificity_tau, 0.99), length.out = 5),
                                               colors = sequential_hcl(n=5, palette = 'Plasma')),
                                           `tau_gt_0.6` = c(`FALSE`='magenta4', `TRUE`='black')
                                  ), show_legend = T, show_annotation_name=F)
    
    H <- H[, rownames(df_pat2nmf)]
    H_viz <- t( scale(H) ) # sample x factor
    pdf(file.path(tmp, 'H.%03d.pdf'), 
        width = 4, height = 7, useDingbats = F, onefile=FALSE)
    draw(Heatmap(t(H), name = 'H', 
                 cluster_columns = FALSE, 
                 cluster_rows = F,
                 right_annotation = row_anno,
                 row_split = pat2nmf[row_i],
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_names_gp = gpar(fontsize = 5),
                 use_raster=T, raster_device = 'CairoPNG')
    )
    draw(Heatmap(H_viz, name = 'H (scaled)', 
                 cluster_columns = FALSE, 
                 right_annotation = row_anno,
                 cluster_rows = F,
                 row_split = pat2nmf[row_i],
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_names_gp = gpar(fontsize = 5),
                 use_raster=T, raster_device = 'CairoPNG')
    )
    dev.off()
    write_rds(pat2nmf, file.path(tmp, 'deliver.patient_best_nmf.rds'))
    p <- ggbarplot(enframe(c(table(pat2nmf))), x='name', y='value', 
                   fill='grey', label=TRUE) +
        labs(y='patient number') + rremove('xlab')
    ggsave(filename = file.path(tmp, 'deliver.patient_best_nmf.pdf'), 
           p, width = 4.5, height = 3, useDingbats = TRUE)

}
#------ Decide the marker genes per factor ------
source('util.nmf.viz.R')
for (r in 4) {
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    W <- read_rds(file.path(tmp, 'W.matrix.rds'))
    # H <- read_rds(file.path(tmp, 'H.matrix.rds'))
    nmf_markers <- define_marker_gene_per_nsnmf(W, n_tolerance = 0)
    str(nmf_markers)
    table(duplicated(unlist(nmf_markers)))
    write_gmx2(nmf_markers, file.path(tmp, 'deliver.nmf_markers.csv'))
    write_rds(nmf_markers, file.path(tmp, 'deliver.nmf_markers.rds'))
    nmf_markers <- lapply(nmf_markers, head, 15)
    dict_nmf2genes <- enframe_list(nmf_markers, name = 'nmf', value = 'gene')
    dict_nmf2genes$gene <- as.character(dict_nmf2genes$gene)
    W_viz <- t(scale( t(W[unique(dict_nmf2genes$gene), ]) ))
    dict_nmf2genes <- dict_nmf2genes[!duplicated(dict_nmf2genes$gene), ]
    W_viz <- W_viz[dict_nmf2genes$gene, ]
    head(rownames(W_viz))
    pdf(file.path(tmp, 'W.%03d.pdf'), 
        width = 4, height = 7.5, useDingbats = F, onefile=F)
    draw(Heatmap(W_viz, name = 'W (scaled)', 
                 cluster_columns = FALSE, 
                 clustering_method_rows = 'complete',
                 cluster_row_slices = TRUE,
                 row_split = deframe(dict_nmf2genes[, c('gene', 'nmf')]), 
                 # row_dend_width  = unit(5, 'mm'),
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_names_gp = gpar(fontsize = 6),
                 use_raster=T, raster_device = 'CairoPNG'))
    draw(Heatmap(W[unique(dict_nmf2genes$gene), ], name = 'W', 
                 cluster_columns = FALSE, 
                 clustering_method_rows = 'complete',
                 cluster_row_slices = TRUE,
                 row_split = deframe(dict_nmf2genes[, c('gene', 'nmf')]), 
                 # row_dend_width  = unit(5, 'mm'),
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_names_gp = gpar(fontsize = 6),
                 use_raster=T, raster_device = 'CairoPNG'))
    
    dev.off()
    
    nmf_markers <- define_marker_gene_per_nsnmf(W, n_tolerance = 2)
    write_gmx2(nmf_markers, file.path(tmp, 'deliver.nmf_markers.tolerant2.csv'))
    write_rds(nmf_markers, file.path(tmp, 'deliver.nmf_markers.tolerant2.rds'))
    
    nmf_markers <- define_marker_gene_per_nsnmf(W, n_tolerance = 1)
    write_gmx2(nmf_markers, file.path(tmp, 'deliver.nmf_markers.tolerant1.csv'))
    write_rds(nmf_markers, file.path(tmp, 'deliver.nmf_markers.tolerant1.rds'))
    
    nmf_top_genes <- top_marker_gene_per_nsnmf(W, topN=30)
    write_gmx2(nmf_top_genes, file.path(tmp, 'deliver.nmf_markers.top.csv'))
    write_rds(nmf_top_genes, file.path(tmp, 'deliver.nmf_markers.top.rds'))
    
    dict_nmf2genes <- enframe_list(nmf_top_genes, 'nmf', 'gene')
    dict_nmf2genes <- dict_nmf2genes[!duplicated(dict_nmf2genes$gene), ]
    W_viz <- t(scale( t(W[unique(dict_nmf2genes$gene), ]) ))
    head(rownames(W_viz))
    pdf(file.path(tmp, 'W.topgene.%03d.pdf'), 
        width = 4, height = 7.5, useDingbats = F, onefile=F)
    draw(Heatmap(W_viz, name = 'W (scaled)', 
                 cluster_columns = FALSE, 
                 clustering_method_rows = 'complete',
                 cluster_row_slices = TRUE,
                 row_split = deframe(dict_nmf2genes[, c('gene', 'nmf')]), 
                 # row_dend_width  = unit(5, 'mm'),
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_names_gp = gpar(fontsize = 3),
                 use_raster=T, raster_device = 'CairoPNG'))
    draw(Heatmap(W[unique(dict_nmf2genes$gene), ], name = 'W', 
                 cluster_columns = FALSE, 
                 clustering_method_rows = 'complete',
                 cluster_row_slices = TRUE,
                 row_split = deframe(dict_nmf2genes[, c('gene', 'nmf')]), 
                 # row_dend_width  = unit(5, 'mm'),
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_names_gp = gpar(fontsize = 3),
                 use_raster=T, raster_device = 'CairoPNG'))
    
    dev.off()

    
}; cat('done.\n')



#------ Survey rank choices ------

# RcppML::mse(A, model$w, model$d, model$h)
reconstruction_err <- c()
for (r in nmf_rank_vec) {
    
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    res <- read_rds( file.path(tmp, 'fastNMF_result_object.rds') )
    err <- RcppML::mse(mat, res$w, res$d, res$h)
    
    reconstruction_err <- c(reconstruction_err, err)
}
df <- data.frame(R=nmf_rank_vec, MSE=reconstruction_err)
write_csv(df, file.path(dir_res, 'survey_ranks.nmf_MSE.csv'))
library(ggplot2)
p <- ggplot( df, aes(x=R, y=MSE )) + geom_point() + geom_line()
ggsave(file.path(dir_res, 'survey_ranks.nmf_MSE.pdf'), p, width = 6, height = 6, useDingbats=F)

cat('done.\n')

#------ Survey rank assignment relations ------
library(clustree)
all_nmf_assignment <- lapply(2:10, function(r){
    o <- read_rds(file.path(dir_res, sprintf('rank%s', r), 'deliver.patient_best_nmf.rds'))
})
names(all_nmf_assignment) <- sprintf('Rank%s', 2:10)
all_obs_names <- names(all_nmf_assignment[[1]])
all_nmf_assignment <- lapply(all_nmf_assignment, function(x) x[all_obs_names])
all_nmf_assignment <- do.call('cbind', all_nmf_assignment)
head(all_nmf_assignment)
p <- clustree(all_nmf_assignment, prefix='Rank')
# print(p)
ggsave(file.path(dir_res, 'survey_ranks.clustree.pdf'), p, width = 7, height = 9, useDingbats=F)

#------ Choose the best K ------
library(factoextra)
library(stringr)
library(ggplot2)
nmf_clustering_finder <- function(k, char_rm) {
    # Given a k, fetch the NMF clustering assignment saved on disk
    # result: a named vector with observation names as names and clustering assignments as values
    o <- file.path('.', 
                   'psbulk', 'fastnmf', sprintf('rank%s', k), 'deliver.patient_best_nmf.rds')
    if (!file.exists(o)) { warning(o, 'not exist!'); return(NA) }
    
    o <- readRDS(o)
    o <- as.numeric( stringr::str_remove_all(o, char_rm) )
    return(o)    
}
nmf_mimic_funcluster <- function(x, k, nmf_clustering_finder, char_rm='fNMF') {
    # x: data matrix (observation x features)                
    # k: number of clusters desired 
    # nmf_clustering_finder: a function fetch the NMF clustering assignment saved on disk
    if (k == 1) {
        clustering <- structure(rep(1, nrow(x)), names=rownames(x))
    } else {
        clustering <- nmf_clustering_finder(k, char_rm=char_rm)
    }
    o <- list()
    o$cluster <- clustering
    return(o)
}

mat <- read_rds(file.path(dir_res, 'input.matrix.rds'))
dim(mat)
for (x in c("silhouette", 'gap_stat')) {    
    pdf(file.path(dir_res, sprintf('survey_ranks.%s.pdf', x)), width = 4,height = 5)
    try(plot(fviz_nbclust(t(mat), 
                          method = x, FUNcluster = nmf_mimic_funcluster, k.max = 10, 
                          nmf_clustering_finder=nmf_clustering_finder, nboot = 50, char_rm='fNMF') +
                 labs(title = 'fNMF', x='rank choice')))
    dev.off()
}
#------ Viz patient-patient overall gene expression to figure out possible clusters ------
mat <- read_rds(file.path(dir_res, 'input.matrix.rds'))
pal_pcr <- c(
    'pCR' = 'seagreen3',
    'RD'='tomato1', 
    'Unknown'='black', 'Excluded'='lightgrey')
idx <- 1:nrow(mat)
pat_cor <- cor(mat[idx, ]); range(pat_cor)
# diag(pat_cor) <- NA
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
col_anno <- read_rds(file.path(dirname(f_in), 'patient_DESeq2.df.rds')); col_anno$patient <- NULL
for (r in nmf_rank_vec) {
    cat(r, '...')
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    pat2nmf <- read_rds( file.path(tmp, 'deliver.patient_best_nmf.rds') )
    col_anno$fNMF <- pat2nmf[rownames(col_anno)]
    fnmf_pal <- structure(rainbow(n=r), names=paste0('fNMF', 1:r))
    pheatmap::pheatmap(
        pat_cor, 
        breaks = c(seq(from=0, to=max(pat_cor[pat_cor!=1]), length.out=9), 1),
        border_color=NA, na_col = NA,
        # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        color = c(colorspace::sequential_hcl(n=10, palette = 'Plasma')),
        cutree_rows=r, cutree_cols=r,
        treeheight_row=15, treeheight_col = 15,
        clustering_method='ward.D2',
        annotation_col  = col_anno, 
        annotation_colors = list(PCR_status=pal_pcr, fNMF=fnmf_pal),
        filename = file.path(tmp, 'pheatmap.cut.sample_correlation.pdf'),
        fontsize = 8, cellwidth=7, cellheight = 7)
    
}
for (x in c("silhouette", 'gap_stat')) {    
    pdf(file.path(dir_res, sprintf('survey_ranks.%s.on_sample_corr.pdf', x)), width = 4,height = 5)
    try(plot(fviz_nbclust(pat_cor, 
                          method = x, FUNcluster = nmf_mimic_funcluster, k.max = 10, 
                          nmf_clustering_finder=nmf_clustering_finder, 
                          nboot = 50, char_rm='fNMF') +
                 labs(title = 'fNMF', x='rank choice')))
    dev.off()
}



#------ Survey ranks similarity ------
source("util.nmf.viz.R")
r <- 2:10    
fin <- file.path(dir_res, sprintf('rank%s', r), 
                 'deliver.nmf_markers.tolerant2.rds')
all(file.exists(fin))
fastnmf_all_ranks <- list()
for (i in seq_along(r)) {
    y <- read_rds(fin[i])
    names(y) <- paste0('Rank', r[i], '_', names(y))
    fastnmf_all_ranks <- c(fastnmf_all_ranks, y)
    rm(y)
}
range(sapply(fastnmf_all_ranks, length))
# fastnmf_all_ranks <- lapply(fastnmf_all_ranks, head, 50)
table(sapply(fastnmf_all_ranks, length)>=10)
gl_R <- fastnmf_all_ranks[sapply(fastnmf_all_ranks, length) >= 10]
name_R <- 'psbulk_rank'
gl_Q <- gl_R
n_genes_R <- sapply(gl_R, length); #n_genes_R <- log10(n_genes_R)
n_genes_Q <- sapply(gl_Q, length); #n_genes_Q <- log10(n_genes_Q)
jaccard_mat <- sapply(gl_R, function(A) {
    sapply(gl_Q, function(B) {
        calc_jaccard(A, B)
    })
})
pval_mat <- sapply(gl_R, function(A) {
    sapply(gl_Q, function(B) {
        calc_enrichment_pval(A, B)
    })
})
qval_mat <- apply(pval_mat, 2, p.adjust, method='bonferroni')
pdf(file.path(dir_res, sprintf('survey_ranks.%s.jaccard.pdf', name_R)), 
    width = 8, height = 6.5, useDingbats = F)
p <- Heatmap(
    jaccard_mat, name='Jaccard', 
    col = heatmap_color_fun_cont(from=0, to=max(jaccard_mat), palette = 'YlOrRd', rev = T), 
    row_dend_width  = unit(5, 'mm'),
    column_dend_height = unit(5, 'mm'),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    clustering_method_rows = 'single',
    clustering_method_columns = 'single',
    # column_split = 4,
    # row_split  = 4,
    height = unit(5, 'inch'),
    width = unit(5, 'inch'),
    right_annotation = rowAnnotation(log10nGene=log10(n_genes_Q)),
    # right_annotation = rowAnnotation(`nGene>=10`=n_genes_Q>=log10(10),
                                     # col=list(`nGene>=10`=c(`TRUE`='black', `FALSE`='lightgrey'))),
    # top_annotation = columnAnnotation(`nGene>=10`=n_genes_R>=log10(10),
    #                                   col=list(`nGene>=10`=c(`TRUE`='black', `FALSE`='lightgrey'))),
    use_raster=T, raster_device = 'CairoPNG')
draw(p)
dev.off()



#------ Visualize genes markers expressions ------
mat <- read_rds(f_expr)
print(dim(mat))
quantile(mat)
rownames(mat)[str_detect(rownames(mat), pattern = 'CLDN')]
cldn_family <- rownames(mat)[str_starts(rownames(mat), pattern = 'CLDN')]

source("util.nmf.viz.R")
pal_pcr <- c(
    'pCR' = 'seagreen3',
    'RD'='tomato1', 
    'Unknown'='black', 'Excluded'='grey')
for (r in 4) {
    cat(r, '...')
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    # nmf_markers <- read_rds(file.path(tmp, 'deliver.nmf_markers.rds'))
    nmf_markers <- read_rds(file.path(tmp, 'deliver.nmf_markers.tolerant2.rds'))
    
    nmf_markers <- lapply(nmf_markers, head, 15)
    dict_nmf2genes <- enframe_list(nmf_markers, name = 'nmf', value = 'gene')
    dict_nmf2genes$gene <- as.character(dict_nmf2genes$gene)
    dict_nmf2genes$nmf <- as.character(dict_nmf2genes$nmf)
    dict_nmf2genes <- dict_nmf2genes[!duplicated(dict_nmf2genes$gene), ]
    pat2nmf <- read_rds( file.path(tmp, 'deliver.patient_best_nmf.rds') )
    table(pat2nmf)
    
    if (r == 5) {
        pat2nmf <- factor(pat2nmf, levels=sprintf('fNMF%s', c(1,4,2,3,5)))
    }
    
    
    col_anno <- read_rds(file.path(dirname(f_in), 'patient_DESeq2.df.rds')); col_anno$patient <- NULL
    col_anno <- col_anno[colnames(mat), ]; dim(col_anno)
    col_anno <- HeatmapAnnotation(
        df=col_anno, 
        col = list(PCR_status=pal_pcr), which = 'column')
    
    manual_genes <- c('EPCAM', 'EGFR', 
                      'KRT5', 'KRT6A', 'KRT6B', 'KRT14', 'KRT17', 'KRT81', 'ACTA2', 'MYLK',
                      'FOXA1', 'KRT19', 
                      'KRT7', 'KRT15', 'KRT16', 'KRT23', 'SLPI', 'LTF', 
                      'KRT8', 'KRT18', 'KRT10', 'AR', 'ESR1', 'PGR')
    manual_genes <- intersect(manual_genes, rownames(mat))
    fibo_genes <- c('LUM', 'FAP',  'MMP3', 'COL1A1', 'COL1A2', 'COL6A1')
    # fibo_genes <- cldn_family
    for (z in c('normed', 'scaled')) {
        row_split <- deframe(dict_nmf2genes[, c('gene', 'nmf')])
        p <- NULL
        p2 <- NULL
        p3 <- NULL
        if (z == 'scaled') {
            mat_viz <- t(scale( t(mat[c(dict_nmf2genes$gene), ]) ))
            mat_viz2 <- t(scale( t( mat[manual_genes, ] )))
            mat_viz3 <- t(scale( t( mat[fibo_genes, ] )))
            col_pal <- heatmap_color_fun_zscore_4
        }
        
        if (z == 'normed') {
            mat_viz <- mat[c(dict_nmf2genes$gene), ]
            mat_viz2 <- mat[manual_genes, ]
            mat_viz3 <- mat[fibo_genes, ]
            
            col_pal <- heatmap_color_fun_cont(
                from=min(c(quantile(mat_viz, .01), quantile(mat_viz2, 0.01), quantile(mat_viz3, 0.01))),
                # from = 0,
                to=max(c(quantile(mat_viz, .99), quantile(mat_viz2, .99), quantile(mat_viz3, .99))))
            # col_pal <- heatmap_color_fun_cont(
            #     from = 0,
            #     to = max(c(max(mat_viz), max(mat_viz2), max(mat_viz3))))
        }

        p3 <- Heatmap(mat_viz3, name = 'fibroblast', 
                      col = col_pal,
                      cluster_columns = F, 
                      cluster_column_slices  = F, 
                      cluster_rows = F,
                      cluster_row_slices = F,
                      column_split = pat2nmf,
                      # row_dend_width  = unit(5, 'mm'),
                      show_row_dend = F, 
                      show_column_dend = F, 
                      column_names_side = 'top',
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 4),
                      border = ifelse(z=='normed', 'white', 'black'), 
                      row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
                      use_raster=T, raster_device = 'CairoPNG')
        
        p2 <- Heatmap(mat_viz2, name = 'markers', 
                      col = col_pal,
                      cluster_columns = F, 
                      cluster_column_slices  = F, 
                      cluster_rows = F,
                      cluster_row_slices = F,
                      column_split = pat2nmf,
                      # row_dend_width  = unit(5, 'mm'),
                      show_row_dend = F, 
                      show_column_dend = F, 
                      column_names_side = 'top',
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 4),
                      border = ifelse(z=='normed', 'white', 'black'), 
                      row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
                      use_raster=T, raster_device = 'CairoPNG')

        p <- Heatmap(mat_viz, name = z, 
                     top_annotation = col_anno,
                     col = col_pal,
                     cluster_columns = F, 
                     cluster_column_slices  = F, 
                     clustering_method_rows = 'ward.D2',
                     clustering_method_columns = 'ward.D2',
                     cluster_row_slices = F,
                     cluster_rows = F,
                     row_split = row_split, 
                     column_split = pat2nmf, 
                     # row_dend_width  = unit(5, 'mm'),
                     show_row_dend = F, 
                     show_column_dend = F, 
                     column_names_side = 'top',
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 4),
                     border = ifelse(z=='normed', 'white', 'black'), 
                     row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
                     use_raster=T, raster_device = 'CairoPNG')
        
        pdf(file.path(tmp, sprintf('expression.%s.pdf', z)), 
            width = 7, height = 9.5, useDingbats = F)
        # draw(p)
        draw(p %v% p2 %v% p3)
        dev.off()
    }
    
}; cat('done\n')

df <- tibble(x=mat['CLDN1', rownames(col_anno)],
       y=structure(col_anno[, 'PCR_status'], names=rownames(col_anno)))
library(ggpubr)
ggboxplot(df, x='y', y='x', add = 'point') + stat_compare_means(comparisons = list(c('pCR', 'RD')))

#------ DEG analysis 1vsRest ------
library(DESeq2)
f_in <- f_expr_deseq2
for (r in 4) {
    cat(r, '...')
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    pat2nmf <- read_rds( file.path(tmp, 'deliver.patient_best_nmf.rds') )
    table(pat2nmf)
    pat2nmf_cat <- gtools::mixedsort( unique(pat2nmf) )
    deg_nmf <- file.path(tmp, 'deg_1vsAll', sprintf('deg_df_%s.rds', pat2nmf_cat))
    names(deg_nmf) <- pat2nmf_cat
    if (!all(file.exists(deg_nmf))) {
        obj <- read_rds(f_in)
        obj$condition <- pat2nmf[as.character(colData(obj)$patient)]
        
        fs::dir_create(file.path(tmp, 'deg_1vsAll'))
        
        for (xx in pat2nmf_cat) {
            cat('DESeq2', xx, '...')
            if (file.exists(file.path(tmp, 'deg_1vsAll', sprintf('deg_df_%s.rds', xx)))) {
                next()
            }
            obj$v <- obj$condition
            idx <- obj$v == xx
            obj$v[idx] <- xx; obj$v[!idx] <- 'rest'
            table(obj$v)
            obj$v <- factor(obj$v)
            obj$v <- relevel(obj$v, ref='rest')
            DESeq2::design(obj) <- formula(~v)
            obj <- DESeq(obj)
            deg <- results(obj)
            
            deg$gene <- rownames(deg)
            deg <- as_tibble(deg)
            deg
            table(deg$padj < 0.01)
            deg <- deg %>% dplyr::arrange(desc(padj < 0.01), desc(log2FoldChange), padj)
            
            # deg_report <- deg %>% dplyr::filter(padj<0.01, log2FoldChange>log2(1))
            deg$fNMF <- xx
            write_csv(deg, file.path(tmp, 'deg_1vsAll', sprintf('deg_df_%s.csv', xx)))
            write_rds(deg, file.path(tmp, 'deg_1vsAll', sprintf('deg_df_%s.rds', xx)))
        }; cat('[done]\n')
    }
    deg_nmf <- lapply(deg_nmf, read_rds)
    deg_nmf_show <- lapply(deg_nmf, function(df) {
        df %>% dplyr::filter(padj < 0.05, log2FoldChange >= log2(1.5)) %>% 
            dplyr::select(gene) %>% deframe()
    })
    write_gmx2(deg_nmf_show, file.path(tmp, 'deg_1vsAll', 'deliver.deg_marker.loose.csv'))
    write_rds(deg_nmf_show, file.path(tmp, 'deg_1vsAll', 'deliver.deg_marker.loose.rds'))
    deg_nmf_show <- lapply(deg_nmf, function(df) {
        df %>% dplyr::filter(padj < 0.01, log2FoldChange >= log2(2)) %>% 
            dplyr::select(gene) %>% deframe()
    })
    write_gmx2(deg_nmf_show, file.path(tmp, 'deg_1vsAll', 'deliver.deg_marker.strigent.csv'))
    write_rds(deg_nmf_show, file.path(tmp, 'deg_1vsAll', 'deliver.deg_marker.strigent.rds'))
    
}

#------ Gene visualization ------
# mat <- read_rds(f_expr_deseq2)
mat <- read_rds(f_expr)
print(dim(mat))
quantile(mat)
source("util.nmf.viz.R")
library(ComplexHeatmap)
pal_pcr <- c(
    'pCR' = 'seagreen3',
    'RD'='tomato1', 
    'Unknown'='black', 'Excluded'='grey')
pal_nmf4 <- c(
    'fNMF1' = '#0066CC',
    'fNMF2' = '#99cc00',
    'fNMF3' = '#ff9933',
    'fNMF4' = '#ff00cc'
); show_col(pal_nmf4)
tapsi_epithelial <- read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/meta/tapsi_epithelial_lineage/epithelial_lineage.rds')
str(tapsi_epithelial)
siyuan_epithelial <- read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/meta/hbca_epithelial_lineage_siyuan/epithelial_lineage.rds')
str(siyuan_epithelial)
y_name <- 'lineage'
hbca_epithelial <- list(
  'Basal' = c('KRT14', 'KRT17', 'DST', 'KRT5', 'SAA1', 'ACTA2', 'SFN',
              'MYLK', 'TAGLN', 'ACTG2'),
  'LumHR' = c('AREG', 'MUCL1', 'AZGP1', 'PIP', 'KRT18', 'AGR2', 'ANKRD30A',
              'S100A14', 'KRT8', 'KRT19'), 
  'LumSec' = c('SCGB2A2', 'SLPI', 'WFDC2', 'LTF', 'KRT15', 'MMP7', 'SCGB3A1', 
               'KRT23', 'CLDN4', 'ALDH1A3')
)

cancer_genes_45 <- read_lines('/volumes/seq/database/CancerGenes/TCGA_breastCancer_genelist_45.txt')
y_name <- 'cancer_gene_45'

keratin <- list(
    Basal=c('KRT5','KRT6A', 'KRT6B', 'KRT14', 'KRT17', 'KRT81'),
    LumSec=c('KRT7', 'KRT15', 'KRT16', 'KRT23'),
    LumHR=c('KRT8', 'KRT18', 'KRT10', 'KRT19'))
y_name <- 'keratin'

misc_clinical <- list(`tnbc`=c('ESR1', 'PGR', 'ERBB2'),
                      `prolif`=c('EPCAM', 'MKI67', 'TOP2A'), 
                      `therapy` = c('EGFR', 'CD274', 'CDK4', 'CDK6'))
y_name <- 'misc_clinical'

cldn_family <- list('Claudins'=cldn_family)
y_name <- 'claudins'
#------ Cell origin ------
tapsi_epithelial <- read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/meta/tapsi_epithelial_lineage/epithelial_lineage.rds')
str(tapsi_epithelial)
cancer_genes_45 <- read_lines('/volumes/seq/database/CancerGenes/TCGA_breastCancer_genelist_45.txt')
y_name <- 'lineage'
y_name <- 'epiLineage'
y_name <- 'misc_clinical'
y_name <- 'hbca_epithelial'
y_name <- 'siyuan_epithelial'
for (r in c(4) ) { ## r=4 is the final solution
    cat(r, '...')
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    pat2nmf <- read_rds( file.path(tmp, 'deliver.patient_best_nmf.rds') )
    table(pat2nmf)
    deg_nmf <- read_rds(file.path(tmp, 'deg_1vsAll', 'deliver.deg_marker.strigent.rds') )
    if (y_name == 'cancer_gene_45') {
        deg_nmf <- read_rds(file.path(tmp, 'deg_1vsAll', 'deliver.deg_marker.loose.rds') )
        print(lapply(deg_nmf, function(x) intersect(cancer_genes_45, x)))
    }
    str(deg_nmf)
    
    if (y_name == c('lineage')) {
        
        gene_viz <- lapply(deg_nmf, function(x) lapply(tapsi_epithelial, function(y) intersect(x, y)))
        
        str(gene_viz)
        gene_viz <- list(LumSec=c(gene_viz$fNMF1$LumSec), 
                         Basal=c(gene_viz$fNMF2$Basal), 
                         LumHR=c(gene_viz$fNMF4$LumHR))
        if (max(sapply(gene_viz, length)) >20){ gene_viz <- lapply(gene_viz, head, 10)}
        gene_viz <- lapply(gene_viz, function(x) intersect(x, rownames(mat)))
        str(gene_viz)
    }    
    if (y_name == 'cancer_gene_45') {
        gene_viz <- list(cg45=cancer_genes_45)
    }
    if (y_name == 'keratin') {
        gene_viz <- keratin
    }
    if (y_name == 'epiLineage') {
        gene_viz <- tapsi_epithelial
    }
    if (y_name == 'misc_clinical') { gene_viz <- misc_clinical }
    if (y_name == 'claudins') { gene_viz <- cldn_family }
    if (y_name == 'hbca_epithelial') { gene_viz <- hbca_epithelial}
    if (y_name == 'siyuan_epithelial') { gene_viz <- siyuan_epithelial; gene_viz <- lapply(gene_viz, head, 15)}
      
    dict_y2genes <- enframe_list(gene_viz, name = 'y', value = 'gene')
    dict_y2genes$gene <- as.character(dict_y2genes$gene)
    dict_y2genes$y <- as.character(dict_y2genes$y)
    dict_y2genes <- dict_y2genes[!duplicated(dict_y2genes$gene), ]
    pat2nmf <- read_rds( file.path(tmp, 'deliver.patient_best_nmf.rds') )
    dict_y2genes <- dict_y2genes[dict_y2genes$gene %in% rownames(mat), ]
    
    
    if ( sum(duplicated(c(unlist(gene_viz)))) != 0 ) {
        xtmp <- c( unlist(gene_viz) )
        xtmp <- xtmp[duplicated(xtmp)]
        gene_viz <- lapply(gene_viz, function(x) x[! x %in% xtmp])
        rm(xtmp)
    }
    str(gene_viz)
        
    
    col_anno <- read_rds(file.path(dirname(f_in), 'patient_DESeq2.df.rds')); col_anno$patient <- NULL
    col_anno <- col_anno[colnames(mat), ]; dim(col_anno)
    col_anno <- HeatmapAnnotation(
        df=col_anno, 
        col = list(PCR_status=pal_pcr), which = 'column')
    
    for (z in c('normed', 'scaled', 'tomax')) {
        row_split <- deframe(dict_y2genes[, c('gene', 'y')])
        p <- NULL
        p2 <- NULL
        p3 <- NULL
        if (z == 'scaled') {
            mat_viz <- t(scale( t(mat[c(dict_y2genes$gene), ]) ))
            col_pal <- heatmap_color_fun_zscore_4
        }
        
        if (z == 'normed') {
            mat_viz <- mat[c(dict_y2genes$gene), ]
            col_pal <- heatmap_color_fun_cont(
                from=quantile(mat_viz, .01),
                to=quantile(mat_viz, .99))
        }
        
        if (z == 'tomax') {
            mat_viz <- t(apply(mat[c(dict_y2genes$gene), ], 1, function(x) x / max(x)))
            col_pal <- heatmap_color_fun_cont(
                from=quantile(mat_viz, .01),
                to=quantile(mat_viz, .99), palette='Lajolla', rev = T)
        }
        
        p <- Heatmap(mat_viz, name = z, 
                     top_annotation = col_anno,
                     col = col_pal,
                     cluster_columns =T,
                     column_split = pat2nmf,
                     cluster_column_slices  = F,
                     
                     clustering_method_columns = 'ward.D2',
                     
                     # cluster_rows = dend_row,
                     cluster_rows = F, 
                     cluster_row_slices = F,
                     clustering_method_rows = 'ward.D2',
                     # row_split = length(unique(row_split)),
                     row_split = row_split,
                     row_dend_width  = unit(5, 'mm'),
                     show_row_dend = T, 
                     show_column_dend = T, 
                     column_names_side = 'top',
                     column_dend_height = unit(5, 'mm'),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 4),
                     border = ifelse(z=='normed', 'white', 'black'), 
                     row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
                     use_raster=T, raster_by_magick = TRUE)
        
        pdf(file.path(tmp, sprintf('expression.%s.%s.pdf', z, y_name)), 
            width = 7, height = 9, useDingbats = F)
        draw(p)
        dev.off()
    }
}
library(ggbeeswarm)

misc_clinical <- list(
`tnbc` = c('ESR1', 'PGR', 'ERBB2'),
`prolif` = c('EPCAM', 'MKI67', 'TOP2A'),
`therapy_gr` = c('AR', 'EGFR', 'FGFR2', 'TACSTD2', 'WEE1'),
`therapy_tme` = c('CD274', 'ANGPT2', 'ANGPT1'),
`therapy_tar1`= c('MYC', 'KRAS', 'AURKA', 'CDK4', 'CDK6'), 
`thearpy_tar3` = c('PIK3CA', 'ATM', 'ATR', 'MDM2', 'MDM4'),
`therapy_tar2` = c('TP53', 'PARP1', 'PARP2', 'BRCA1', 'BRCA2')
)
theme_set(theme_pubr())
for (r in c(4) ) { ## r=4 is the final solution
  cat(r, '...')
  tmp <- file.path(dir_res, sprintf('rank%s', r))
  pat2nmf <- read_rds( file.path(tmp, 'deliver.patient_best_nmf.rds') )
  pat2nmf <- enframe(pat2nmf, name = 'patient', value = 'NMF4')
  unlist(misc_clinical)
  misc_clinical_gnames <- c("ESR1", "PGR", "ERBB2", "AR", "EPCAM", "MKI67", "TOP2A", "EGFR" , "CD274", "CDK4" , "CDK6" )
  # for (i in seq_along(misc_clinical)) {
    # misc_clinical_gnames <- misc_clinical[[i]]
  misc_clinical_gnames <- unlist(misc_clinical)
  misc_clinical_gnames[which(!misc_clinical_gnames %in% rownames(mat))]
    # boxplot(lapply(misc_clinical_gnames, function(x) mat[x, ]))
    df <- mat[as.character(misc_clinical_gnames), ] %>% as.data.frame() %>%
      rownames_to_column('gname') %>%
      gather(., 'patient', 'value', colnames(mat), factor_key=F)
    
    col_anno <- read_rds(file.path(dirname(f_in), 'patient_DESeq2.df.rds')); col_anno$patient <- NULL
    col_anno <- col_anno[colnames(mat), ]
    col_anno <- rownames_to_column(col_anno, 'patient')
    df <- left_join(df, col_anno, by='patient')
    df <- left_join(df, pat2nmf, by='patient')
    head(df)
        
    library(ggbeeswarm)
    df$gname <- factor(df$gname, levels = misc_clinical_gnames) ## important to avoid bug
    
    library(ggpubr); library(rstatix)
    
    p <- ggboxplot(df, x='gname', y='value', color = 'NMF4', palette = pal_nmf4, 
                   outlier.shape = NA, fill=NA) + 
      geom_quasirandom(aes(color=NMF4), size=0.4, dodge.width=.8) +
      rotate_x_text(90) +
      labs(x='', y='expression (log2CPM)')

    p
    stat.test <- df %>%
      group_by(gname) %>%
      rstatix::wilcox_test(value ~ NMF4, p.adjust.method='fdr') #%>%
      # adjust_pvalue(method = "BH") %>%
      # add_significance("p.adj")
    head(stat.test)
    stat.test$p.adj.txt <- signif(stat.test$p.adj, digits = 3)
    View(stat.test)
    stat.test <- stat.test %>%
      add_xy_position(x = "gname", dodge = 0.8)
    p <- p + stat_pvalue_manual(
      stat.test,  label = 'p.adj.txt', tip.length = 0.007, 
      bracket.nudge.y = 0.5)
    ggsave(filename = file.path(tmp, sprintf('expression.boxplot.clinical_genes_%s.pdf', 'combo')),
           plot =  p, width = 0.8 * length(misc_clinical_gnames),
           height = 5)
  # }
    
}


#------ Enrichment on marker genes only ------
source('pkg_enricher.R')
fs::dir_create(dir_res)

for (r in 4) {
    cat(r, '... ')
    tmp <- file.path(dir_res, sprintf('rank%s', r))
    
    # list_module_val <- read_rds(file.path(dir_res_module, 'module_content.list.rds'))
    # str(list_module_val)
    list_module_val <- read_rds(file.path(tmp, 'deliver.nmf_markers.rds'))
    list_module_val <- read_rds(file.path(tmp, 'deliver.nmf_markers.tolerant2.rds'))
    str(list_module_val)
    list_module_val <- lapply(list_module_val, function(G){
        G <- gsub(pattern = 'MT-', replacement = 'MT', x=G) # MT-ATP8: MTATP8
        G <- gsub(pattern = "\\.\\d+", replacement='', x=G) # POLR2J3.1: POLR2J3
        G <- unlist(lapply(G, function(g) {
            g_matched <- limma::alias2Symbol(g, species = "Hs")
            if (purrr::is_empty(g_matched)) { return(g) }
            return(g_matched)
        }))
        G <- .gene_symbol2entrezID(x = G, db=org.Hs.eg.db)
        return(G)
    })
    str(list_module_val)
    module_names <- names(list_module_val); str(module_names)
    #------ MsigDb Hallmark ------
    
    for (d in c('H', 'C2', 'C5')) {
        db_use <- switch (d,
                          H = db_msigdb_H,
                          C2 = db_msigdb_C2,
                          C5 = db_msigdb_C5
        ) 
        db_name <- switch (d,
                           H = 'MSigDB Hallmark',
                           C2 = 'MSigDB C2 Curated',
                           C5 = 'MSigDB C5 GO'
        )
        message(db_name)
        if (!file.exists(file.path(tmp, sprintf('enricher.%s.rds', db_name)))) {
            df_enricher <- read_rds(file.path(tmp, sprintf('enricher.%s.rds', db_name)))
        } else {
            df_enricher <- suppressMessages(run_enricher(
                X = list_module_val, y=db_use,
                pvalueCutoff = 0.05,
                pAdjustMethod = 'BH', maxGSSize=5000, verbose = T))
            write_csv(
                df_enricher, 
                file.path(tmp, sprintf('enricher.%s.csv', db_name)))
            write_rds(
                df_enricher, 
                file.path(tmp, sprintf('enricher.%s.rds', db_name)))
        }    
        df_enricher <- dplyr::filter(df_enricher, p.adjust < 0.05)
        df_enricher$genesymbol <- sapply(df_enricher$geneID, function(x) {
            x <- unlist(str_split(x, '/'))
            x <- .gene_entrezID2symbol(x, db=org.Hs.eg.db)
            return( paste(x, collapse = '/'))
        })
        df_enricher$Enrichment <- -log10(df_enricher$p.adjust)
        if (nrow(df_enricher) == 0) {next()}
        print(table(df_enricher$query_name))
        df_enricher <- df_enricher %>% 
            dplyr::arrange(query_name, desc(Enrichment))
        write_csv(
            df_enricher, 
            file.path(tmp, sprintf('enricher.%s.report.csv', db_name)))
        unique(df_enricher$query_name)

        df_enricher_report <- lapply(levels(df_enricher$query_name), function(x) {
            y <- df_enricher[df_enricher$query_name == x, ]
            y <- y[order(y$Enrichment, decreasing = T), ]
            return(y$ID)
        }); names(df_enricher_report) <- levels(df_enricher$query_name); str(df_enricher_report)
        write_gmx2(df_enricher_report, file.path(tmp, sprintf('enricher.%s.gmx.csv', db_name)))
        df_enricher <- df_enricher %>%
            dplyr::group_by(query_name) %>%
            dplyr::top_n(n=15, wt=Enrichment) %>% dplyr::ungroup()
        db_sgntr_names <- sort(unique(df_enricher[, 'ID', drop=T]))
        mat_enricher <- tidy_df2mat(
            df=df_enricher, row_str = 'ID', col_str = 'query_name',
            val_str = 'Enrichment',
            row_levels = db_sgntr_names, col_levels = module_names,
            fill_na = 0)
        rownames(mat_enricher) <- gsub(pattern = 'HALLMARK_', replacement = '', rownames(mat_enricher))
        print(dim(mat_enricher))
        
        pdf(file.path(tmp, 
                      sprintf('enricher.heatmap.anno_modules.%s.pdf', db_name)),
            width = 3+4, height = ifelse(d=='H', 5, 15))
        draw(heatmap_pval(mat_enricher = mat_enricher, row_title=db_name,
                          width = unit(3, 'inch'),qval_cutoff = 0.05,
                          row_names_gp = gpar(fontsize = 6)))
        dev.off()
        
    }
    
}