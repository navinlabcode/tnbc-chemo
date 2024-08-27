## Helper functions related to gene enrichment
## Yun Yan


suppressPackageStartupMessages({
    # library(future, lib.loc = "/usr/lib64/R/library")
    library(future)
    library(future.apply)
    # library(future.apply, lib.loc = "/usr/lib64/R/library")
  library(tidyverse); library(ggplot2); library(ggpubr); library(patchwork)
  library(cli); library(tictoc); library(glue); library(scales); library(tools)
  library(ComplexHeatmap); library(circlize); library(dendextend); library(dendsort)
    library(clusterProfiler) # 3.14.2
    library(org.Hs.eg.db)
  # library(org.Mm.eg.db)
  library(msigdbr)

  .gene_symbol2entrezID <- function(x, db){
    # x <- c('Hopx', 'Sftpc', 'Sftpb', 'yunyan')
    # .gene_symbol2entrezID(x=x, db=org.Mm.eg.db)
    #    Hopx   Sftpc   Sftpb  yunyan
    # "74318" "20389" "20388"      NA
    mapIds(x=db, keys=x, column='ENTREZID', keytype='SYMBOL')
  }
  .gene_entrezID2symbol <- function(x, db){
    mapIds(x=db, keys=x, column='SYMBOL', keytype='ENTREZID')
  }
  gmt2list <- function(x){
    res <- read_lines(x, skip_empty_rows = T)
    res_names <- sapply(res, function(s){
      as.character(unlist(stringr::str_split(s, pattern ='\t', n=2)))[1]
    })
    res <- lapply(res, function(s){
      s <- as.character(unlist(stringr::str_split(s, pattern ='\t')))
      s <- stringr::str_trim(s)
      sort(s[3:length(s)]) ## exclude 1st and 2nd slot
    })
    names(res) <- res_names
    return(res)
  }
  list_enframe <- function(L){
    L_size <- sapply(L, length)
    res <- data.frame(
      name=rep(names(L), L_size), 
      value=unlist(L), stringsAsFactors = F)
    res$name <- factor(as.character(res$name), levels=names(L))
    return(res)
  }
})

db_msigdb_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
db_msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
db_msigdb_C5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory='BP') %>% 
  dplyr::select(gs_name, entrez_gene)
db_msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'CP') %>% 
    dplyr::select(gs_name, entrez_gene)

pathwaysH = split(x = db_msigdb_H$entrez_gene, f = db_msigdb_H$gs_name)
pathwaysC2 = split(x = db_msigdb_C2$entrez_gene, f = db_msigdb_C2$gs_name)
pathwaysC5 = split(x = db_msigdb_C5$entrez_gene, f = db_msigdb_C5$gs_name)
# module_entrezid_list <- lapply(module_genes_list, .gene_symbol2entrezID, db=org.Hs.eg.db)
# print(sapply(module_entrezid_list, function(v) v[is.na(v)]))
# module_entrezid_list <- lapply(module_entrezid_list, function(z) unique(as.numeric(z[!is.na(z)])))


library(future.apply)
#' Wrapper to run clusterProfiler::enricher for multiple batches of genes.
#' 
#' @param y A data frame. See the parameter `TERM2GENE` in the function `clusterProfiler::enricher`.
#' @param ... Options to run `clusterProfiler::enricher`
#' @param X A list of entrez IDs.
#' @param verbose 
#' @param use_future_apply 
#'
#' @return A data frame of `enrichResult` with an additional column 'query_name'.
#' @author Yun Yan
#' @examples 
#' run_enricher(
#'   list('gene_set_A'=c(123, 34), 'gene_set_B'=c(100, 300)),
#'   MsigDB_H)
run_enricher <- function(X, y, verbose=FALSE, use_future_apply=TRUE, ...){
  x_names <- names(X)
  # func <- lapply
  # if (use_future_apply) {func <- future.apply::future_lapply}
  df_enrich <- future.apply::future_lapply(x_names, function(x){
    if (verbose) {cat(x, '... ')}
    o <- clusterProfiler::enricher(gene=X[[x]], TERM2GENE=y, ...)
    # minGSSize=20, maxGSSize=2000,
    # pAdjustMethod='BH', pvalueCutoff = 0.05
    if (purrr::is_empty(o)) {return(NULL)}
    o <- o@result
    if (nrow(o)==0) {return(NULL)}
    o$query_name <- x
    return(o)
  }); if (verbose) {cat('[done]\n')}
  df_enrich <- do.call('rbind', df_enrich)
  df_enrich$query_name <- factor(df_enrich$query_name, levels=x_names)
  return(df_enrich)
}


#' Wrapper to run clusterProfiler::GSEA for multiple batches of genes.
#'
#' @param X A list of named numeric vector where names are the EntrezID and values are logFC.
#' @param y A data frame. See the parameter `TERM2GENE` in the function `clusterProfiler::GSEA`.
#' @param verbose 
#' @param use_future_apply 
#' @param ... 
#'
#' @return A data frame of `enrichResult` with an additional column 'query_name'.
#' @export
#' @author Yun Yan
#' @examples
#' data(geneList, package="DOSE")
run_gsea <- function(X, y, verbose=FALSE, use_future_apply=F, ...){
  x_names <- names(X)
  func <- lapply
  if (use_future_apply) {func <- future.apply::future_lapply}
  df_gsea <- func(x_names, function(x){
    if (verbose) {cat(x, '... ')}
    o <- clusterProfiler::GSEA(geneList=X[[x]], TERM2GENE=y, ...)
    # minGSSize=20, maxGSSize=2000,
    # pAdjustMethod='BH', pvalueCutoff = 0.05
    if (purrr::is_empty(o)) {return(NULL)}
    o <- o@result
    if (nrow(o)==0) {return(NULL)}
    o$query_name <- x
    return(o)
  }); if (verbose) {cat('[done]\n')}
  df_gsea <- do.call('rbind', df_gsea)
  df_gsea$query_name <- factor(df_gsea$query_name, levels=x_names)
  return(df_gsea)
}


#------ 2. Long df => Matrix ------
tidy_df2mat <- function(
  df, row_str, col_str, val_str,
  row_levels, col_levels, fill_na=0){
  
  i <- match(as.character(df[, row_str, drop=T]), row_levels)
  j <- match(as.character(df[, col_str, drop=T]), col_levels)
  v <- df[, val_str, drop=T]
  mat <- Matrix::sparseMatrix(
    i=i, j=j, x=v, 
    dims=c(length(row_levels), length(col_levels)))
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- fill_na
  rownames(mat) <- row_levels
  colnames(mat) <- col_levels
  return(mat)
}
# "ID" "p.adjust" "Count" "GeneRatio" "query_name" 
#------ 3. P.adjust Heatmap ------
# "ID" "p.adjust" "Count" "GeneRatio" "query_name" eval(parse(text='1/10'))
heatmap_pval <- function(mat_enricher, qval_cutoff=0.01, ...){
  p <- Heatmap(
    mat_enricher, name='-log10(qval)',
    cluster_rows = T, cluster_columns = F,
    show_row_names = T, show_column_names = T,
    clustering_method_rows  = 'complete',
    clustering_method_columns = 'complete',
    show_row_dend = F,
    row_names_side = 'left', column_names_side = 'top',
    row_dend_width = unit(5, 'mm'),
    column_dend_height = unit(5, 'mm'),
    column_names_rot = 45,
    
    col = colorRamp2(c(0, -log10(qval_cutoff), max(mat_enricher)),
                     c('black', 'white', 'firebrick')),
    # width = unit(5, 'inch'),
    ...
  )
  return(p)
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

deframe_to_fgsea_genes <- function(df, z=1, y=2, x=3){
    # z: which cluster
    # y: names
    # x: values
    
    stopifnot(ncol(df) == 3)
    df <- as.data.frame(df)
    lv <- gtools::mixedsort(unique(df[, z]))
    res <- lapply(lv, function(lvx) {
        structure(df[df[, z] == lvx, x], names=df[df[, z] == lvx, y])
    })
    names(res) <- lv
    return(res)
}