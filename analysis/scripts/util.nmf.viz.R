#
# Helper functions related to NMF 
# - visualization
# - processing
#
# Author: Yun Yan
#


library(ComplexHeatmap); library(circlize); library(colorspace)
library(dendextend)
library(dendsort)
library(Matrix)
library(scales)

scale_nmf_gene_loading <- function(X){
    # X: genes x factors
    res <- apply(X, 2, function(v) scales::rescale(x=v, to=c(0,1)))
    rownames(res) <- rownames(X)
    colnames(res) <- colnames(X)
    return(res)
}
zscore_nmf_gene_loading <- function(X){
    # X: genes x factors
    res <- apply(X, 2, function(v) scale(v))
    rownames(res) <- rownames(X)
    colnames(res) <- colnames(X)
    return(res)
}

define_marker_gene_per_nsnmf <- function(X, n_tolerance=0){
    # X: NMF gene loading matrix (W) genes x factors
    # Ref: Reuben and Itai Yanai's ST paper.
    ## Don't use rank=1 because of ties. 
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
    factor_names <- colnames(X)
    res <- lapply(factor_names, function(J){
        v <- sort(X[, J], decreasing = TRUE)
        gnames <- names(v)
        return( head(gnames, topN) )
    })
    names(res) <- factor_names
    return(res)
}


define_marker_gene_per_nmf_by_ttest <- function(x, INDEX, do.pval.adj=TRUE, significance_thres=0.05){
    # x: gene loading matrix (genes, factors)
    # INDEX: named vector. factor:membership
    levs <- gtools::mixedsort(unique(INDEX))
    rnames <- rownames(x)
    cnames <- colnames(x)
    
    calc_single_t_test_pval <- function(a, b){
        obj <- t.test(x=a, y=b, alternative = 'greater')
        return(obj$p.value)
    }
    
    res <- lapply(levs, function(l){
        j_focus <- INDEX == l
        j_other <- INDEX != l
        pvals <- apply(x, 1, function(v){
            calc_single_t_test_pval(v[j_focus], v[j_other])
        })
        if (do.pval.adj) {
            pvals <- p.adjust(pvals, method = 'fdr')
        }
        is_good <- pvals < significance_thres
        return(rnames[is_good])
    })
    names(res) <- levs
    return(res)
}

heatmap_color_num_fun <- colorRamp2(
    seq(from=-1, to=1, length.out = 5),
    colors = diverge_hcl(n=5, palette = 'Blue-Red 3')
)
heatmap_color_fun <- heatmap_color_num_fun

heatmap_cor_mat <- function(
    cor_mat, hclust_obj, hclust_color_k=4,
    run_name='X', topdf='heatmap_corrmat.pdf', cor_name='Corr', 
    cluster_rows_cols_slice = F, ...) {
    
    cluster_rows_cols <- F
    # cluster_rows_cols_slice <- F
    if (!is.null(hclust_obj)){
        hclust_obj <-  dendextend::color_branches(
            hclust_obj, k = hclust_color_k, col = ruok::init_random_color(1:hclust_color_k))
        cluster_rows_cols <- hclust_obj
        cluster_rows_cols_slice <- TRUE
    }
    
    p1 <- Heatmap(
        cor_mat, name=cor_name, 
        # cluster_rows = F,
        # cluster_columns = F,
        cluster_rows = cluster_rows_cols,
        cluster_columns = cluster_rows_cols,
        cluster_row_slices = cluster_rows_cols_slice,
        cluster_column_slices = cluster_rows_cols_slice,
        row_dend_width  = unit(5, 'mm'),
        column_dend_height = unit(5, 'mm'),
        
        show_row_names = F, 
        show_column_names = F,
        
        # row_title = sprintf('%d rows', nrow(cor_mat)),
        # column_title  = sprintf('%d columns', ncol(cor_mat)),
        
        height = unit(5, 'inch'),
        width = unit(5, 'inch'),
        na_col = 'cyan', ...
    )
    
    pdf(topdf, height = 6.5, width = 6.5)
    draw(p1,
         column_title=sprintf('%s (%d columns)', run_name, ncol(cor_mat)),
         column_title_side='bottom')
    dev.off()
    return(p1)
}
heatmap_color_correlation_fun <- colorRamp2(
    seq(from=-1, to=1, length.out = 5),
    colors = diverge_hcl(n=5, palette = 'Blue-Red 3')
)

heatmap_color_fun_zscore_4 <- colorRamp2(
    seq(from=-4, to=4, length.out = 5),
    colors = diverge_hcl(n=5, palette = 'Blue-Red 3')
)
heatmap_color_fun_zscore_2 <- colorRamp2(
    seq(from=-2, to=2, length.out = 5),
    colors = diverge_hcl(n=5, palette = 'Blue-Red 3')
)


heatmap_color_fun_zscore <- heatmap_color_fun_zscore_4

heatmap_color_fun_cont <- function(from, to, palette='Viridis', rev=F) {
    colorRamp2(
        seq(from=from, to=to, length.out = 5),
        colors = sequential_hcl(n=5, palette = palette, rev = rev))
}

get_cell_order <- function(obj, y){
    structure( tibble::deframe(FetchData(obj, y)), names=Cells(obj) )
}

write_gmx2 <- function (L, path) 
{
    column_names <- names(L)
    if (purrr::is_empty(column_names)) {
        column_names <- as.character(seq_along(L))
    }
    n_row <- max(sapply(L, length))
    n_col <- length(L)
    df <- as.data.frame(matrix(NA, nrow = n_row, ncol = n_col))
    for (j in seq_len(n_col)) {
        i_to <- length(L[[j]])
        if (i_to == 0) { next() }
        df[1:i_to, j] <- L[[j]]
    }
    colnames(df) <- column_names
    readr::write_csv(x = df, path = path, col_names = TRUE, na = "")
}


calc_jaccard <- function(a, b) {
    length(unique(intersect(a, b))) / length(unique(union(a, b)))
}
calc_enrichment_pval <- function(a, b, N=24705) {
    # phyper(A&B-1, B, ALL-B, A, lower.tail= FALSE)
    phyper(length(intersect(a, b)) - 1, 
           length(b), 
           N - length(b), 
           length(a), 
           lower.tail = FALSE)
}

self_pairwise_run <- function(x, func, run_parallel = TRUE, ...) {
    n_rows <- n_cols <- length(x)
    x_names <- names(x)
    
    mat_ij <- lapply(1:n_rows, function(i){
        Jx <- i:n_rows
        Ix <- rep(i, length(Jx))
        return( cbind(i=Ix, j=Jx) )
    })
    
    mat_ij <- do.call(rbind, mat_ij)
    if (run_parallel) {
        v <- future_apply(mat_ij, MARGIN = 1, FUN = function(r){
            func(x[[r[1]]], x[[r[2]]], ...)
        })
    } else {
        v <- apply(mat_ij, MARGIN = 1, FUN = function(r){
            func(x[[r[1]]], x[[r[2]]], ...)
        })
    }
    mat_ij <- cbind(mat_ij, v=v)
    
    mat_ij <- mat_ij[mat_ij[, 'v'] != 0, ] ## only keep non-zero slots -> sparse matrix    
    
    mat_ij_lower_tri <- cbind(i = mat_ij[, 'j'], j=mat_ij[, 'i'], v=mat_ij[, 'v'])
    mat_ij_lower_tri <- mat_ij_lower_tri[mat_ij_lower_tri[, 'i'] != mat_ij_lower_tri[, 'j'], ]
    
    mat_ij <- rbind(mat_ij, mat_ij_lower_tri)
    # print(mat_ij)
    mat_ij <- new('dgTMatrix', 
                  i=as.integer(mat_ij[, 'i']-1),
                  j=as.integer(mat_ij[, 'j']-1),
                  x=as.numeric(mat_ij[, 'v']),
                  Dim = c(n_rows, n_cols))
    
    rownames(mat_ij) <- colnames(mat_ij) <- x_names
    return(mat_ij)
}

