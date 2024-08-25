#!/usr/bin/Rscript
#--------------------------
# Perform fastNMF with a fixed rank
# 
# Yun Yan (yun.yan@uth.tmc.edu)
#--------------------------
suppressPackageStartupMessages({
    library(RcppML)
    library(Biobase)
    library(matrixStats)
    library(tictoc); library(cli)
    library(readr); library(tidyverse)
})

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 0) {
    f_in         <- args[1]
    nmf_rank     <- args[2]
    dir_res      <- args[3]
    nmf_rank <- as.numeric(nmf_rank)
    stopifnot( length(nmf_rank) == 1 )
} else {
    ## example
    f_in <- './expr_matrix/each_sample/ARTC44.gc.NonNegCenterMat.rds'
    nmf_rank <- 1e3
    dir_res <- './metamodule_fnmf/nmf_each_sample/ARTC44/K14'
}
fs::dir_create(dir_res)
nmf_rank <- as.numeric(nmf_rank)

# reconstruction_err.txt <= this is a must because it will recognized by the caller script
# W.matrix.rds
# H.matrix.rds
 ## 
#------ Preprare nmf input matrix ------

mat <- read_rds(f_in)
mat <- as.matrix(mat); print(class(mat))
cat(sprintf('Input matrix = %s', f_in), 'with dimension: \n')
print(dim(mat))
cat(sprintf('Rank = %s', nmf_rank), '. \n')
cat(sprintf('Out dir = %s', dir_res), '. \n')
fs::dir_create(dir_res)


for (r in nmf_rank) {
    if (r > ncol(mat)) {
        cat(sprintf('Ranks is greater than number of observations. %s ranks v.s. %s observations.', r, ncol(mat)))
        cat('[skipped]\n')
        err <- -1024
        write_lines(err, file.path(dir_res, 'reconstruction_err.txt'))
        quit(save='no')
    }
}
nmf_rank <- nmf_rank[nmf_rank <= ncol(mat)]

#------ Run NMF ------
for (r in nmf_rank) {
    
    if (file.exists(file.path(dir_res, 'fastNMF_result_object.rds'))) {
        cat('load existing fastNMF object...\n')
        res <- read_rds(file.path(dir_res, 'fastNMF_result_object.rds'))
    } else {
        res <- nmf(A=mat, k=r, tol=1e-5, L1 = c(0.01, 0.01), seed = 42)
        write_rds( res, file.path(dir_res, 'fastNMF_result_object.rds') )
    }
    print(dim(res$w)) # genes x factor
    print(dim(res$h)) # factor x samples
    W <- res$w
    H <- res$h
    # W[1:3, 1:3]
    # H[1:3, 1:3]
    
    rownames(W) <- rownames(mat)
    colnames(H) <- colnames(mat)
    colnames(W) <- rownames(H) <- paste0('fNMF', 1:r)
    
    write_rds(W, file.path(dir_res, 'W.matrix.rds'))
    write_rds(H, file.path(dir_res, 'H.matrix.rds'))
    
}

suppressPackageStartupMessages({
    library(ruok)
    library(ComplexHeatmap)
    library(ggpubr)
    library(ggplot2)
})
#------ Decide the dictionary sample ~ NMF ------
for (r in nmf_rank) {
    H <- read_rds(file.path(dir_res, 'H.matrix.rds')) # factors x sample
    pat2nmf <- rownames(H)[apply(H, 2, which.max)]
    names(pat2nmf) <- colnames(H)
    H_viz <- t( scale(H) ) # samples x factors
    pdf(file.path(dir_res, sprintf('H.pdf')), 
        width = 4, height = 7, useDingbats = F)
    draw(Heatmap(H_viz, name = 'H', 
                 cluster_columns = FALSE, 
                 clustering_method_rows = 'complete',
                 cluster_row_slices = TRUE,
                 row_split = pat2nmf,
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 show_row_names = F,
                 row_title_rot = 0, column_title_rot=90,
                 row_names_gp = gpar(fontsize = 5),
                 use_raster=T, raster_device = 'CairoPNG')
    )
    dev.off()
    write_rds(pat2nmf, file.path(dir_res, 'deliver.cell_best_nmf.rds'))
    p <- ggbarplot(enframe(c(table(pat2nmf))), x='name', y='value', 
                   fill='grey', label=TRUE) +
        labs(y='cell number') + rremove('xlab')
    ggsave(filename = file.path(dir_res, 'deliver.cell_best_nmf.pdf'), 
           p, width = 4.5, height = 3, useDingbats = TRUE)
    
    
}
#------ Decide the marker genes per factor ------
define_marker_gene_per_nsnmf <- function(X, n_tolerance=2){
    # X: NMF gene loading matrix (W) genes x factors
    # Ref: Reuben and Itai Yanai's ST paper.
    ## Don't naively use rank=1 because some have ties. 
    factor_names <- colnames(X)
    mat_row_ismax <- t( apply(X, 1, function(xx) xx == max(xx, na.rm=TRUE)) )
    colnames(mat_row_ismax) <- factor_names
    res <- lapply(factor_names, function(J){
        n_violate <- 0
        v <- sort(X[, J], decreasing = TRUE, na.last=TRUE)
        gnames <- names(v)
        ismax <- mat_row_ismax[gnames, J]
        out <- c()
        i <- 1
        while( i <= length(gnames) ){
            if (!ismax[i]) { n_violate <- n_violate + 1}
            if (n_violate > n_tolerance) { break() }
            if (ismax[i] ) {out <- c(out, gnames[i])}
            i <- i + 1
        }
        return(out)
    })
    names(res) <- factor_names
    return(res)
}

top_marker_gene_per_nsnmf <- function(X, topN=50) {
    # X: NMF gene loading matrix (W) genes x factors
    # Ref: This follows the Reuben and Itai Yanai's ST paper.
    factor_names <- colnames(X)
    res <- lapply(factor_names, function(J){
        v <- sort(X[, J], decreasing = TRUE)
        gnames <- names(v)
        return( head(gnames, topN) )
    })
    names(res) <- factor_names
    return(res)
}

for (r in nmf_rank) {
    W <- read_rds(file.path(dir_res, 'W.matrix.rds'))
    # H <- read_rds(file.path(dir_res, 'H.matrix.rds'))
    nmf_markers <- define_marker_gene_per_nsnmf(W, n_tolerance = 0)
    nmf_markers <- nmf_markers[sapply(nmf_markers, function(xx) length(xx) > 0)]
    str(nmf_markers)
    table(duplicated(unlist(nmf_markers)))
    write_gmx2(nmf_markers, file.path(dir_res, 'deliver.nmf_markers.csv'))
    write_rds(nmf_markers, file.path(dir_res, 'deliver.nmf_markers.rds'))
    nmf_markers <- lapply(nmf_markers, head, 15)
    dict_nmf2genes <- enframe_list(nmf_markers, name = 'nmf', value = 'gene')
    dict_nmf2genes$gene <- as.character(dict_nmf2genes$gene)
    W_viz <- t(scale( t(W[unique(dict_nmf2genes$gene), ]) ))
    dict_nmf2genes <- dict_nmf2genes[!duplicated(dict_nmf2genes$gene), ]
    W_viz <- W_viz[dict_nmf2genes$gene, ]
    head(rownames(W_viz))
    pdf(file.path(dir_res, sprintf('W.pdf')), 
        width = 7, height = 10, useDingbats = F)
    draw(Heatmap(W_viz, name = 'W', 
                 cluster_columns = FALSE, 
                 clustering_method_rows = 'complete',
                 cluster_row_slices = TRUE,
                 row_split = deframe(dict_nmf2genes[, c('gene', 'nmf')]), 
                 # row_dend_width  = unit(5, 'mm'),
                 show_row_dend = F, 
                 show_column_dend = F, 
                 column_names_side = 'top',
                 row_title_rot = 0, column_title_rot=90,
                 row_names_gp = gpar(fontsize = 5),
                 use_raster=T, raster_device = 'CairoPNG'))
    dev.off()
    
    nmf_markers <- define_marker_gene_per_nsnmf(W, n_tolerance = 2)
    write_gmx2(nmf_markers, file.path(dir_res, 'deliver.nmf_markers.tolerant2.csv'))
    write_rds(nmf_markers, file.path(dir_res, 'deliver.nmf_markers.tolerant2.rds'))
    
    nmf_top_genes <- top_marker_gene_per_nsnmf(W, topN=50)
    write_gmx2(nmf_top_genes, file.path(dir_res, 'deliver.nmf_markers.top.csv'))
    write_rds(nmf_top_genes, file.path(dir_res, 'deliver.nmf_markers.top.rds'))
    
}; cat('done.\n')

#------ Visualize genes markers expressions ------
print(dim(mat))

for (r in nmf_rank) {
    cat(r, '...')
    nmf_markers <- read_rds(file.path(dir_res, 'deliver.nmf_markers.rds'))
    nmf_markers <- lapply(nmf_markers, head, 15)
    dict_nmf2genes <- enframe_list(nmf_markers, name = 'nmf', value = 'gene')
    dict_nmf2genes$gene <- as.character(dict_nmf2genes$gene)
    dict_nmf2genes$nmf <- as.character(dict_nmf2genes$nmf)
    dict_nmf2genes <- dict_nmf2genes[!duplicated(dict_nmf2genes$gene), ]

    for (z in c('normed', 'scaled')) {
        row_split <- deframe(dict_nmf2genes[, c('gene', 'nmf')])
        p <- NULL
        p2 <- NULL
        if (z == 'scaled') {
            mat_viz <- t(scale( t(mat[dict_nmf2genes$gene, ]) ))
            mat_viz2 <- t(scale( t( mat[manual_genes, ] )))
            col_pal <- heatmap_color_fun_zscore_4
        }
        
        if (z == 'normed') {
            manual_genes <- c('EPCAM', 
                              'KRT5', 'KRT6A', 'KRT17', 'KRT23', 'KRT81', 
                              'TP63', 'MME',
                              'KRT7', 'KRT8', 'KRT18', 'KRT19', 'FOXA1'
            )
            manual_genes <- intersect(manual_genes, rownames(mat))
            
            mat_viz <- mat[c(dict_nmf2genes$gene), ]
            mat_viz2 <- mat[manual_genes, ]
            
            col_pal <- heatmap_color_fun_cont(
                from=min(c(quantile(mat_viz, .05), quantile(mat_viz2, .05))),
                to=max(c(quantile(mat_viz, .95), quantile(mat_viz2, .95))))
        }
        
        p2 <- Heatmap(mat_viz2, name = 'manual', 
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
                      row_names_gp = gpar(fontsize = 4),
                      column_names_gp = gpar(fontsize = 4),
                      show_column_names = F,
                      row_title_rot = 0, column_title_rot=90,
                      border = T, 
                      row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
                      use_raster=T, raster_device = 'CairoPNG')
        
        p <- Heatmap(mat_viz, name = z, 
                     # top_annotation = col_anno,
                     col = col_pal,
                     cluster_columns = F, 
                     cluster_column_slices  = T, 
                     clustering_method_rows = 'ward.D2',
                     clustering_method_columns = 'complete',
                     cluster_row_slices = F,
                     row_split = row_split, 
                     column_split = pat2nmf, 
                     # row_dend_width  = unit(5, 'mm'),
                     show_row_dend = F, 
                     show_column_dend = F, 
                     column_names_side = 'top',
                     row_names_gp = gpar(fontsize = 4),
                     column_names_gp = gpar(fontsize = 4),
                     show_column_names = F,
                     row_title_rot = 0, column_title_rot=90,
                     border = T, 
                     row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
                     use_raster=T, raster_device = 'CairoPNG')
        
        pdf(file.path(dir_res, sprintf('expression.%s.pdf', z)), 
            width = 7, height = 8.5, useDingbats = F)
        # draw(p)
        draw(p %v% p2)
        dev.off()
    }
    
}; cat('done\n')

#------ Survey rank choices ------
reconstruction_err <- c()
for (r in nmf_rank) {
    res <- read_rds( file.path(dir_res, 'fastNMF_result_object.rds') )
    err <- RcppML::mse(mat, res$w, res$d, res$h)
    
    write_lines(err, file.path(dir_res, 'reconstruction_err.txt'))
    reconstruction_err <- c(reconstruction_err, err)
}


cat('done.\n')