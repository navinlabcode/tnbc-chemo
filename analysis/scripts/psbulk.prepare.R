## Prepare psbulk expression data
##
## Yun Yan

library(Seurat)
library(Matrix)
library(tidyverse)
library(fs)
library(DESeq2)
dir_proj <- './'
args <- commandArgs(trailingOnly = TRUE)

#-------------------------- INPUT --------------------------  

#------ Tumor only compartment ------
if (!is_empty(args)) {
    fpath_sr3_ready <- args[[1]]
} else {
    fpath_sr3_ready <- file.path(dir_proj, 'ready.sr3.rds')
}

message(fpath_sr3_ready)

#------ Outdir ------
stopifnot(file.exists(fpath_sr3_ready))
dir_res <- file.path(dir_proj, 'psbulk'); fs::dir_create(dir_res)

#-------------------------- START --------------------------  

sr3 <- read_rds(fpath_sr3_ready)
print(sr3)
length(unique(as.character(sr3$patient))) # 99 -> 97 patients

if ( ! file.exists( file.path(dir_res, 'normalized.DESeq2VST.matrix.rds') )) {
    # ps-bulk RNA by patients
    expr_mat <- GetAssayData(sr3, assay = 'RNA', slot = 'counts') %>% as.matrix()
    expr_mat <- mat_tapply(t(expr_mat), as.character(sr3$patient), sum) # sample x genes
    expr_mat <- t(expr_mat)
    # print(expr_mat[1:3, 1:3])
    write_rds(expr_mat, file.path(dir_res, 'counts.matrix.rds'))
    # expr_mat <- read_rds( file.path(dir_res, 'counts.matrix.rds') )
    print(dim(expr_mat))
    
    # df <- unique(sr3[[]][, c('patient', 'PCR_status', 'RCB_status')])
    df <- unique(sr3[[]][, c('patient', 'pCR_status', 'RCB_status', 'archetype')])
    rownames(df) <- df$patient
    df <- df[colnames(expr_mat), ]
    dds <- DESeqDataSetFromMatrix(
        countData = expr_mat, 
        colData = df, 
        design = ~1 # do not specifying a model
    )
    dds_norm <- vst(dds)
    vst_mat <- assay(dds_norm)
    # identical(dds_norm$patient, names(deframe(enframe(table(sr3$patient)))))
    dds_norm$log10nCells <- log10(c(table(sr3$patient))[dds_norm$patient])
    write_rds(dds, file.path(dir_res, 'DESeq2_obj.rds'))
    write_rds(dds_norm, file.path(dir_res, 'DESeq2_normed_obj.rds'))
    write_rds(vst_mat, file.path(dir_res, 'normalized.DESeq2VST.matrix.rds'))
    write_rds(as.data.frame(colData(dds_norm)), file.path(dir_res, 'patient_DESeq2.df.rds'))
    expr_mat <- t( expr_mat )
    lib_size <- rowSums(expr_mat)
    expr_mat <- expr_mat / lib_size * 1e6
    expr_mat <- log2(expr_mat + 1)
    expr_mat <- t( expr_mat )
    write_rds(expr_mat, file.path(dir_res, 'normalized.matrix.rds'))
    write_rds(df, file.path(dir_res, 'patient_metainfo.rds'))
}
cat('DONE')


#------ Viz patient-patient correlation ------
    
# dir_res <- file.path('.//psbulk')
mat <- read_rds( file.path(dir_res, 'normalized.DESeq2VST.matrix.rds') )
# mat <- t( scale(t(mat)) )
print(dim(mat))

g_var <- rowVars(mat)
g_avg <- rowMeans(mat)

plot(g_avg, g_var)
idx <- 1:nrow(mat)
pat_cor <- cor(mat[idx, ]); range(pat_cor)
write_rds(pat_cor, file.path(dir_res, 'pat_cor_matrix.rds'))
library(ComplexHeatmap)
library(pheatmap)
identical(rownames(pat_cor), rownames(colData(dds_norm)))
annotation_row <- as.data.frame(colData(dds_norm))
annotation_row$patient <- NULL
pheatmap::pheatmap(
    pat_cor, 
    treeheight_row=15, treeheight_col = 15,
    clustering_method='ward.D',
    annotation_row = annotation_row, 
    fontsize = 8, cellwidth=10, cellheight = 10,
    filename = file.path(dir_res, 'pheatmap.no_cut.sample_correlation.pdf'))
for (cut_k in c(2, 3, 4, 5, 6, 7, 8)) {
    pheatmap::pheatmap(
        pat_cor, 
        cutree_rows=cut_k, cutree_cols=cut_k,
        annotation_row = annotation_row, 
        treeheight_row=15, treeheight_col = 15,
        clustering_method='ward.D',
        fontsize = 8, 
        cellwidth=10, cellheight = 10,
        filename = file.path(dir_res, sprintf('pheatmap.cut_%s.sample_correlation.pdf', cut_k)))
}

