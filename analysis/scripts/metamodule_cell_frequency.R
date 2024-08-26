## Module score and cellular frequency of metaprograms
## - AddModuleScore
## - Cell frequency
## - Comparison of pCR and RD
##      - module score
##      - cell freq
## - Comparison of the four archetypes
##      - module score
##      - cell freq
##
## Yun Yan (yun.yan@uth.tmc.edu)
##

suppressPackageStartupMessages({
    library(Seurat); library(S4Vectors)
    library(tidyverse); library(ggplot2); library(ggpubr); library(patchwork)
    theme_set(theme_pubr(base_size = 12, legend='right')); library(ruok)
    library(cli); library(tictoc); library(glue); library(scales); library(tools)
    library(ComplexHeatmap); library(circlize)
    library(purrr); library(stringr); library(readr)
    library(msigdbr)
    my_scatter_theme <- theme_pubr(base_size = 12, legend = 'right') %+replace%
        theme(axis.text=element_text(size=rel(.2)), axis.title=element_text(size=rel(.3)), axis.ticks=element_line(size = rel(.3)), panel.border = element_rect(fill = NA, size=rel(.5)), axis.line = element_blank())
    my_scatter_themevoid <- theme_pubr(base_size = 12, legend = 'bottom') %+replace% theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), panel.border = element_blank(), axis.line = element_blank())  
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
    source('util.nmf.viz.R')
})


#-------------------------- Inputs --------------------------    
dir_res <- '.'
fpath_sr3_ready <- file.path(dir_res, 'ready.sr3.rds')

f_cellmeta_tumor <- ''
dr_use <- 'UMAP' # 'tSNE' / 'UMAP'
snn_res_viz <- NULL

sr3 <- read_rds(fpath_sr3_ready) 
df_cellmeta <- read_rds(file.path(dir_res, 'sr3_metadata.df.rds'))

#-------------------------- Gene signatures --------------------------    
if (TRUE) {
    f_mm_memb <- file.path(dir_res, 'metamodule_fnmf', 'all_ranks', 
                           'jaccard', 'manual_alt_clean_wardD2', 
                           'jaccard-manual.metamodule_membership.csv')
    module_content <- read_rds(
        file.path('.',
                  'metamodule_fnmf', 
                  'MM_alt_clean_byscore_wardD2', 
                  'deliver.mm_markers.rds')
    )
    str(module_content)
    module_focuse <- setdiff(names(module_content), 'M14'); cat(module_focuse)
    module_focuse_bio <- c('G2/M', 'Mitochondria', 'Ribosome', 'Stress', 'Interferon',
                           'HLA','S/G1', 'Hypoxia', 'Basal','EMT',
                           'LumSec', 'Cholestero', 'StressER')
    names(module_focuse_bio) <- module_focuse
    publish_top_markers <- 20
    module_content_use <- lapply(module_content, function(x) head(x, publish_top_markers))
    module_content_use <- module_content_use[module_focuse]
    str(module_content_use)
    names(module_content_use) <- paste0(
        names(module_content_use), '__',
        make.names(module_focuse_bio[names(module_content_use)]))
    # db <- c(module_content_use, db_msigdb_H)
    db <- module_content_use
    dir_res_viz <- file.path(
      dir_res, 
      'viz_signature_MM_alt_clean_byscore_wardD2'); fs::dir_create(dir_res_viz)
}

## Other gene signatures
## for example: immunogenic: cGAS-sting and CIN
if (FALSE) {
  dir_res_viz <- file.path(dir_res, 'viz_signature_immunogenic'); fs::dir_create(dir_res_viz)
  ## c-gas/sting (STING1=TMEM173)
  db_cgas_sting <- read_lines('/volumes/USR1/yyan/project/tnbc_pre_atlas/lib/signature_cGAS_STING.txt')
  ## IFN alpha + beta
  db_ifn_A <- read_lines('/volumes/USR1/yyan/project/tnbc_pre_atlas/lib/signature_IFN_alpha.txt')
  db_ifn_B <- read_lines('/volumes/USR1/yyan/project/tnbc_pre_atlas/lib/signature_IFN_beta.txt')
  db_ifn1 <- c(db_ifn_A, db_ifn_B)
  ## IFN gamma
  db_ifn2 <- read_lines('/volumes/USR1/yyan/project/tnbc_pre_atlas/lib/signature_IFN_gamma.txt')
  ## Genome instability CIN25
  db_cin25 <- read_lines('/volumes/USR1/yyan/project/tnbc_pre_atlas/lib/signature_CIN25.txt')
  db <- list(CGAS=db_cgas_sting, IFN1=db_ifn1, IFN2=db_ifn2, CIN25 = db_cin25)
}


str(db)
unlist(db)[which(!unlist(db) %in% rownames(sr3))]
#-------------------------- Calculating --------------------------    

db_name_prefix <- 'Signature_'
if (! is.null(db_name_prefix) ) {
    names(db) <- paste0(db_name_prefix, names(db))
}
names(db) <- paste0(names(db), '_Seurat')
if ( file.exists(file.path(dir_res_viz, 'addmodulescore.df.rds')) ) {
    df_signature_res <- read_rds(file.path(dir_res_viz, 'addmodulescore.df.rds'))
    stopifnot(identical(rownames(df_signature_res), Cells(sr3)))
    sr3 <- AddMetaData(sr3, metadata = df_signature_res)
    
} else {
    sr3 <- AddModuleScore(sr3, features = db, name = names(db))
    for (i in seq_along(db)) {
        x <- names(db)[i]
        y <- sprintf('%s%s', x, i)
        cat(sprintf('%s', x), '\n')
        sr3[[x]] <- sr3[[y]]
        sr3[[y]] <- NULL
    }
    colnames(sr3[[]])
    # df_signature_res <- sr3[[]][, grepl(db_name_prefix, colnames(sr3[[]])) ]
    df_signature_res <- sr3[[]][, names(db) ]
    print(colnames(df_signature_res))
    write.csv(x=df_signature_res, file = file.path(dir_res_viz, 'addmodulescore.df.csv'))
    write_rds(x=df_signature_res, file.path(dir_res_viz, 'addmodulescore.df.rds'))
}
cat('Done AddModuleScore.\n')
#stop('quit')
# q(save='no')
#-------------------------- viz --------------------------    

pal_pcr_old <- c(
    'pCR' = 'seagreen3',
    'RD'='tomato1', 
    'Unknown'='black', 'Excluded'='lightgrey')
pal_pcr <- c('pCR'='#53D43F',
             'RD'='#811C9A', 
             'Unknown'='black', 
             'Excluded'='grey',
             'Removed'='ghostwhite')
pal_nmf4 <- c(
    'ARC1' = '#0066CC',
    'ARC2' = '#99cc00',
    'ARC3' = '#ff9933',
    'ARC4' = '#ff00cc'
)

#------ z: patient groupgs info ------
patient_info_pcr  <- sr3[[]][, c('patient', 'PCR_status')] %>% unique() %>% deframe()
patient_info_nmf4 <- read_rds(file.path(dir_res, "psbulk/fastnmf/rank4/deliver.patient_best_nmf.rds"))

sr3$PCR_label <- sr3$PCR_status
cells_in_use <- sr3$PCR_status %in% c('pCR', 'RD')

sr3$NMF4 <- patient_info_nmf4[as.character(sr3$patient)]

#-------------------------- Viz each patient --------------------------
mm_signature_names <- colnames(df_signature_res)
print(mm_signature_names)

fs::dir_create(file.path(dir_res_viz, 'viz_each_sample'))
pat_avail <- unique(sr3$patient); str(pat_avail)

## Only shows those patients with MM-calculated

mm_memb <- read_csv(f_mm_memb)
pat_avail <- unique(mm_memb$sample); str(pat_avail)

# for (run_i in pat_avail) {
for (run_i in c('ARTC100', 'ARTC95')) {
    cat(run_i, '...')
    # df_signature_res
    cells_in_use <- Cells(sr3)[sr3$patient == run_i]
    if (length(cells_in_use) > 1e3 | length(cells_in_use) < 100) { next() } ## don't visualize big/small data
    df_signature_use <- df_signature_res[cells_in_use, ]
    
    df_signature_use <- df_signature_use[, mm_signature_names]
    dim(df_signature_use)
    
    ## score, zscore, binary, binary_v2 matrices are to be calcualted and visualized
    mat_run_i <- as.matrix(df_signature_use) # cells x features
    mat_run_zscore_i <- scale(mat_run_i)
    mat_run_binary_i <- mat_run_i >= 0.1
    
    mat_run_i <- t(mat_run_i)
    mat_run_zscore_i <- t(mat_run_zscore_i)  # features x cells
    mat_run_binary_i <- t(mat_run_binary_i) * 1
    
    mat_run_binary_i_v2 <- mat_run_binary_i * mat_run_zscore_i
    mat_run_binary_i_v2 <- mat_run_binary_i_v2 > 0
    mat_run_binary_i_v2 <- mat_run_binary_i_v2 * 1
    
    stopifnot(identical(rownames(mat_run_i), rownames(mat_run_zscore_i)))
    mm_signature_names_pretty <- str_remove_all(rownames(mat_run_i), 'Signature_') %>% str_remove_all(., '_Seurat')
    rownames(mat_run_i) <- rownames(mat_run_zscore_i) <- rownames(mat_run_binary_i) <- rownames(mat_run_binary_i_v2) <- mm_signature_names_pretty
    
    ## ordering cells based on max zscore -- for visualization purpose
    cells_split <- apply(mat_run_zscore_i, 2, function(o) {
        if (max(o) <= 0) {return('Unresolved')}
        names(o)[nnet::which.is.max(o)]
    }); print(table(cells_split)); print(table(cells_split)/length(cells_split))
    cell_idx <- order(cells_split)
    mat_run_i <- mat_run_i[, cell_idx]
    mat_run_zscore_i <- mat_run_zscore_i[, cell_idx]
    mat_run_binary_i <- mat_run_binary_i[, cell_idx]
    mat_run_binary_i_v2 <- mat_run_binary_i_v2[, cell_idx]
    cells_split <- cells_split[cell_idx]
    
    pal_mm <- structure(
        c(ArchR::paletteDiscrete(mm_signature_names_pretty), 'Grey'),
        names=c(mm_signature_names_pretty, 'Unresolved')) #; show_col(pval_mm)
    ht_cell_anno <- columnAnnotation(
        meta_module = cells_split, 
        col = list(meta_module=pal_mm), 
        annotation_legend_param = list(ncol=3, legend_direction='horizontal', grid_height=unit(2, 'mm')),
        show_legend = T, show_annotation_name=F
    )
    head(rownames(mat_run_i),5)
    heatmap_color_fun_zscore <- colorRamp2(
        seq(from=-4, to=4, length.out = 5),
        colors = diverge_hcl(n=5, palette = 'Blue-Red 3')
    )
    heatmap_color_fun_score <- colorRamp2(
        seq(from=quantile(mat_run_i, 0.01), to=quantile(mat_run_i, 0.99), length.out = 5),
        colors = sequential_hcl(n=5, palette = 'Plasma')
    )
    heatmap_color_fun_binary <- c(`1`='black', `0`='lightgrey')
    heatmap_color_fun_binary_v2 <- c(`1`='darkblue', `0`='lightgrey')
    for (mat_type in c('score', 'zscore', 'binary', 'binary_v2')) {
        if (mat_type == 'score' ) {tmp_mat <- mat_run_i; tmp_col <- heatmap_color_fun_score}
        if (mat_type == 'zscore' ) {tmp_mat <- mat_run_zscore_i; tmp_col <- heatmap_color_fun_zscore}
        if (mat_type == 'binary' ) {tmp_mat <- mat_run_binary_i; tmp_col <- heatmap_color_fun_binary}
        if (mat_type == 'binary_v2' ) {tmp_mat <- mat_run_binary_i_v2; tmp_col <- heatmap_color_fun_binary_v2}
        ht_obj <- Heatmap(
            tmp_mat, name=mat_type, 
            col=tmp_col,
            column_split = cells_split,
            
            cluster_rows = F, 
            cluster_row_slices = T, 
            show_row_dend = F, show_column_dend = F,
            
            cluster_columns = F, 
            
            show_column_names = F,
            show_row_names = T,
            
            row_title_rot = 0, 
            column_title_rot=0,
            row_gap = unit(0, 'mm'), 
            column_gap  = unit(0, 'mm'), 
            border = T,
            top_annotation = ht_cell_anno,
            row_names_gp = gpar(fontsize = 6),
            # left_annotation = ht_gene_anno2,
            
            column_title = sprintf("%s (%s cells)", run_i, 
                                   scales::comma(ncol(tmp_mat))), 
            use_raster=TRUE, 
            raster_by_magick = TRUE,
            heatmap_height = unit(2, 'inch')
            # annotation_legend_param = list(
            #     direction = "horizontal")
        )
        pdf(file.path(dir_res_viz, 'viz_each_sample', 
                      sprintf('demo_MM.%s.%s.pdf', run_i, mat_type)), 
            width = 5, height = 3)
        draw(ht_obj, heatmap_legend_side = "right",
             annotation_legend_side = 'bottom')
        dev.off()
    }
    cat('done.')
}; cat('Done\n')


#-------------------------- Gene Score diff --------------------------    

#------ violin ------
library(ggbeeswarm)
for (z in c('nmf4', 'pcr')) {
    if (z == 'pcr') {grouping_z <- 'PCR_label'; pal_use <- pal_pcr; p_width <- 2}
    if (z == 'nmf4') {grouping_z <- 'NMF4'; pal_use <- pal_nmf4; p_width <- 2/2*4}
    sr3$grouping_z <- sr3[[grouping_z]]
    for (x in colnames(df_signature_res)) {
        p <- VlnPlot(sr3,
                     features = x, 
                     group.by='grouping_z', pt.size = 0) +
            labs(title=str_remove(x, db_name_prefix) %>% str_remove('HALLMARK_'))
        tmp_xlab <- sr3[[]] %>% dplyr::group_by(grouping_z) %>% 
            dplyr::summarise_at(x, mean) %>% deframe() %>%
            formatC(digits = 2) %>% pretty_table2str()
        p <- p + 
            stat_compare_means() +
            rremove('xlab') + rremove('legend') + 
            labs(y='signature score') +
            rotate_x_text(0) 
        p <- p + scale_x_discrete(labels=tmp_xlab)
        if (z == 'pcr'){p <- p + scale_x_discrete(limits=c('pCR', 'RD'), labels=tmp_xlab[c('pCR', 'RD')]) }
        p <- p + scale_fill_manual(values=pal_use) +
          theme(axis.ticks.length=unit(0.07,"inch"))
        ggsave(file.path(dir_res_viz,
                         sprintf('VlnPlot.by%s.%s.pdf', z, x) ),
               p, width = p_width, height = 3)
    }
    p_list <- list()
    for (x in colnames(df_signature_res)) {    
        summary_bypat <- sr3[[]] %>% 
            dplyr::group_by(grouping_z, patient) %>% 
            dplyr::summarise_at(x, mean)
        p <- ggplot(summary_bypat, aes_string(x='grouping_z', y=x)) +
          # geom_beeswarm(aes(color=grouping_z), cex=2) + 
          geom_boxplot(outlier.shape=NA) +
          geom_quasirandom(aes(color=grouping_z), size=2, alpha=0.5, pch=16) + 
          rremove('y.title') +
          labs(title=str_remove(x, db_name_prefix) %>% str_remove('HALLMARK_'),
               y='signature score')
        if (z == 'pcr'){p <- p + scale_x_discrete(limits=c('pCR', 'RD')) }
        p <- p + 
          ggpubr::stat_mean(pch=5, cex=2, color="gold") +
          stat_compare_means(aes(label = paste0("p=", after_stat(p.format)))) +
            rremove('xlab') + rremove('legend') +
            rotate_x_text(0) +
            scale_color_manual(values=pal_use) +
          theme_pubr(base_size = 6) +
          theme(axis.ticks.length=unit(0.07,"inch")) + 
          rremove('legend') + rremove('x.title')
        p_list <- c(p_list, list(p))
        ggsave(file.path(dir_res_viz,
                         sprintf('VlnPlot.by%s.%s.patient.pdf', z, x) ),
               p, width = p_width, height = 3)
    }
    # p_list[[1]]
    ggsave(file.path(dir_res_viz,
                     sprintf('VlnPlot.by%s.%s.patient.pdf', z, 'compact') ),
           ggarrange(plotlist = p_list, nrow=1), width = p_width * length(p_list), height = 3, limitsize=F, useDingbats = F)
    
}


#------ heamtap ------
    
mat_signature_score <- sr3[[]][, grepl('_Seurat', colnames(sr3[[]]) )] %>% 
    as.matrix()
dim(mat_signature_score) # cells x features
mat_signature_score_bypat <- mat_tapply(
    mat = mat_signature_score, INDEX = as.character(sr3$patient), 
    FUN = function(xx) {mean(xx, na.rm=T)})
dim(mat_signature_score_bypat)
write_rds(mat_signature_score_bypat, file.path(dir_res_viz, 'mat_patient_modulescore.rds'))
colnames(mat_signature_score_bypat)
## Only show the patients used for finding meta-modules -- those patients have sufficient tumor cells
mat_viz <- t(mat_signature_score_bypat[, mm_signature_names])
mat_viz <- mat_viz[, pat_avail]
dim(mat_viz)
# Heatmap(mat_viz)
col_anno <- HeatmapAnnotation(
    df= data.frame(
        NMF4=patient_info_nmf4[colnames(mat_viz)],
        PCR_status=patient_info_pcr[colnames(mat_viz)]),
    col = list(PCR_status=pal_pcr, NMF4=pal_nmf4), which = 'column')

for (ii in c(0, 1, 2, 3, 4)) {
    ## 0: clustering and coloring are independent
    ## 1: coloring prioritize pCR/RD info then NMF4. Clustering is done per slice.
    ## 2: coloring prioritize NMF4 then pCR/RD.  Clustering is done per slice.
    ## 3: coloring prioritize pCR/RD only.
    ## 4: coloring prioritize NMF4 only.
    print(ii)
    col_split <- NULL
    if (ii == 0) {col_split <- NULL} 
    if (ii == 1) {col_split <- paste0(
        patient_info_pcr[colnames(mat_viz)],
        patient_info_nmf4[colnames(mat_viz)], sep='_')}
    if (ii == 2) {col_split <- paste0(
        patient_info_nmf4[colnames(mat_viz)], 
        patient_info_pcr[colnames(mat_viz)], sep='_')}
    if (ii == 3) {col_split <- patient_info_pcr[colnames(mat_viz)]}
    if (ii == 4) {col_split <- patient_info_nmf4[colnames(mat_viz)]}
    heatmap_color_fun_score <- colorRamp2(
        seq(from=quantile(mat_viz, 0.01), to=quantile(mat_viz, 0.99), length.out = 5),
        colors = sequential_hcl(n=5, palette = 'Plasma')
    )
    p <- Heatmap(
        mat_viz, name = 'score',
        col = heatmap_color_fun_score,
        show_row_names = TRUE,
        show_column_names = TRUE,
        top_annotation = col_anno,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        column_split = col_split,
        
        cluster_columns = ifelse(ii == 0, TRUE, FALSE), 
        clustering_method_columns = 'ward.D2',
        cluster_column_slices = T,
        column_title = NULL,
        border = T,
        row_dend_width  = unit(5, 'mm'),
        column_dend_height  = unit(5, 'mm'),
        row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
        # use_raster=T, raster_device = 'CairoPNG'
    )
    pdf(file.path(dir_res_viz, sprintf('heatmap.signature_score.%s.pdf', ii)), 
        width = 7.5, height = 4, useDingbats = F)
    draw(p)
    dev.off()
}

#-------------------------- Cell fraction diff --------------------------    

pal_pcr2_old <- c(`pCR`='seagreen3', `RD`='tomato1', `n.s.`='lightgrey')
pal_pcr2 <- c(`pCR`='#53D43F', `RD`='#811C9A', `n.s.`='lightgrey')
pat_avail <- unique(mm_memb$sample); str(pat_avail)


cells_in_use <- c(sr3$patient %in% pat_avail); print(table(cells_in_use))

df_cellmeta_use <- sr3[[]][cells_in_use, ]
mat_signature_score <- sr3[[]][cells_in_use, grepl('_Seurat', colnames(sr3[[]]) )] %>% 
    as.matrix()
dim(mat_signature_score) # cells x features
mat_signature_score <- mat_signature_score[, mm_signature_names] ## ad-hoc: only check YY's modules
colnames(mat_signature_score) <- str_remove(colnames(mat_signature_score), 'Signature_') %>% str_remove_all(., '_Seurat')
print(colnames(mat_signature_score))

classifier_method <- c( 'hybrid', 'winner_v2', 'hybrid_v2', 'hybrid_v3', 'hybrid_v4')
library(ComplexHeatmap)
mat_signature_zscore <- scale(mat_signature_score) # cells x features
source("util.nmf.viz.R")
library(ggbeeswarm)
for (cfm in classifier_method) {
    fs::dir_create(file.path(dir_res_viz, cfm))
    cat(cfm, '\n')
    
    #------ Decide assignment for each cell ------
    ## [in use] hybrid: a cell can use multiple mm
    if (cfm == 'hybrid') {        
        mat_signature_binary <- (mat_signature_score >= 0.1)
        mat_signature_binary <- mat_signature_binary * 1
        unresolved_c <- apply(mat_signature_binary, 1, function(v) all(v==0))
        print(table(unresolved_c))
        unresolved_v <- rep(0, length(unresolved_c))
        unresolved_v[unresolved_c] <- 1
        mat_signature_binary <- cbind(mat_signature_binary, Unresolved=unresolved_v)
    }
    ## [DEPRECATED] winner: a cell only use the highest mm
    if (cfm == 'winner') { 
        ## Problem: What if a MM is overall higher than other MMs, then this MM is always the winner, e.g., the ribo
        mat_signature_binary <- apply(mat_signature_score, 1, function(v) {
            v == max(v) 
        }) %>% t()
    }
    ## [DEPRECATED] winner: a cell only use the highest zscored mm
    if (cfm == 'winner_v2') {
        # exclude mito and ribo MM
        mm_use <- setdiff(colnames(mat_signature_score), 
                          c('M02__Mitochondria', 'M03__Ribosome'))        
        pat_avail <- unique(df_cellmeta_use$patient); str(pat_avail)
        mat_signature_binary <- list()
        for (run_i in pat_avail) {
            cells_i <- rownames(df_cellmeta_use)[df_cellmeta_use$patient == run_i]
            mat_i <- mat_signature_score[cells_i, mm_use, drop=F]
            mat_i <- scale(mat_i)
            mat_i_binary <- apply(mat_i, 1, function(v) {
                res <- rep(0, length(v))
                ## A cell has no meta-modules -- Unresolved
                if (max(v) <= 0) {return( rep(NA, length(v)) )} 
                res[nnet::which.is.max(v)] <- 1
                return(res)
            })
            mat_i_binary <- t(mat_i_binary)
            colnames(mat_i_binary) <- colnames(mat_i)
            unresolved_c <- apply(mat_i_binary, 1, function(v) all(is.na(v)))
            unresolved_v <- rep(0, length(unresolved_c))
            unresolved_v[unresolved_c] <- 1
            mat_i_binary <- cbind(mat_i_binary, Unresolved=unresolved_v)
            mat_i_binary[is.na(mat_i_binary)] <- 0
            mat_signature_binary <- c(mat_signature_binary, list(mat_i_binary))
            
        }
        mat_signature_binary <- do.call('rbind', mat_signature_binary)
        mat_signature_binary <- mat_signature_binary[rownames(mat_signature_score), ]

    }
    ## [DEPRECATED] hybrid_v2: a cell can use multiple mm as long as the zcored MM >=0.5 
    if (cfm == 'hybrid_v2') {
        cfm_hybrid_cutoff <- 0.5
        mm_use <- colnames(mat_signature_score)
        pat_avail <- unique(df_cellmeta_use$patient); str(pat_avail)
        mat_signature_binary <- list()
        for (run_i in pat_avail) {
            cells_i <- rownames(df_cellmeta_use)[df_cellmeta_use$patient == run_i]
            mat_i <- mat_signature_score[cells_i, mm_use, drop=F]
            mat_i <- scale(mat_i)
            mat_i_binary <- apply(mat_i, 1, function(v) {
                res <- rep(0, length(v))
                ## A cell has no meta-modules -- Unresolved
                if (max(v) <= 0) {return( rep(NA, length(v)) )} 
                res <- c(v >= cfm_hybrid_cutoff) * 1
                return(res)
            })
            mat_i_binary <- t(mat_i_binary)
            colnames(mat_i_binary) <- colnames(mat_i)
            unresolved_c <- apply(mat_i_binary, 1, function(v) all(is.na(v)))
            unresolved_v <- rep(0, length(unresolved_c))
            unresolved_v[unresolved_c] <- 1
            mat_i_binary <- cbind(mat_i_binary, Unresolved=unresolved_v)
            mat_i_binary[is.na(mat_i_binary)] <- 0
            mat_signature_binary <- c(mat_signature_binary, list(mat_i_binary))
            
        }
        mat_signature_binary <- do.call('rbind', mat_signature_binary)
        mat_signature_binary <- mat_signature_binary[rownames(mat_signature_score), ]
        
        
    }
    ## [DEPRECATED] hybrid_v3: a cell can use multiple mm as long as the zcored MM > 0 and score >= 0.1
    if (cfm == 'hybrid_v3') {
        cfm_hybrid_score_cutoff <- 0.1
        cfm_hybrid_zscore_threshold <- 0
        mm_use <- colnames(mat_signature_score)
        pat_avail <- unique(df_cellmeta_use$patient); str(pat_avail)
        mat_signature_binary <- list()
        for (run_i in pat_avail) {
            cells_i <- rownames(df_cellmeta_use)[df_cellmeta_use$patient == run_i]
            mat_i <- mat_signature_score[cells_i, mm_use, drop=F]
            mat_zscore_i <- scale(mat_i)
            
            mat_i_binary <- (mat_i >= cfm_hybrid_score_cutoff) * (mat_zscore_i > cfm_hybrid_zscore_threshold)
            mat_i_binary <- mat_i_binary * 1
            
            colnames(mat_i_binary) <- colnames(mat_i)
            unresolved_c <- apply(mat_i_binary, 1, function(v) all(v==0))
            unresolved_v <- rep(0, length(unresolved_c))
            unresolved_v[unresolved_c] <- 1
            mat_i_binary <- cbind(mat_i_binary, Unresolved=unresolved_v)
            mat_i_binary[is.na(mat_i_binary)] <- 0
            mat_signature_binary <- c(mat_signature_binary, list(mat_i_binary))
            
        }
        mat_signature_binary <- do.call('rbind', mat_signature_binary)
        mat_signature_binary <- mat_signature_binary[rownames(mat_signature_score), ]
        
        
    }
    ## [DEPRECATED] hybrid_v4: a cell can use multiople mm as long as the glocal zcored MM > 0
    if (cfm == 'hybrid_v4') {
        cfm_hybrid_score_cutoff <- 0.1
        cfm_hybrid_zscore_threshold <- 0
        mat_signature_binary <- (mat_signature_zscore > cfm_hybrid_zscore_threshold) * 1
        mat_signature_binary <- mat_signature_binary[rownames(mat_signature_zscore), ]
    
    }
    
    #------ Visualization ------
    mat_patient_cellpct <- mat_tapply(
        mat = mat_signature_binary, INDEX = df_cellmeta_use$patient, 
        FUN = function(xx) {sum(xx, na.rm = TRUE) / length(xx)})
    # sort(rowSums(mat_patient_cellpct))
    colnames(mat_patient_cellpct)
    # rownames(mat_patient_cellpct) <- gtools::mixedsort(unique(df_cellmeta_use$patient)) # patients x features
    df_patient_cellpct <- as.data.frame(mat_patient_cellpct)
    df_patient_cellpct$patient <- rownames(mat_patient_cellpct)
    # df: signature - patient - pct - grouping
    df_patient_cellpct <- tidyr::gather(
        df_patient_cellpct, 'signature', 'pct', 
        colnames(mat_patient_cellpct))
    head(df_patient_cellpct)
    
    df_patient_cellpct$PCR <- patient_info_pcr[as.character(df_patient_cellpct$patient)]
    df_patient_cellpct$NMF4 <- patient_info_nmf4[as.character(df_patient_cellpct$patient)]
    write_rds(df_patient_cellpct, 
              file.path(dir_res_viz, cfm, 'df_patient_cellpct.rds'))
    write_rds(mat_patient_cellpct, 
              file.path(dir_res_viz, cfm, 'mat_patient_cellpct.rds'))
    write_rds(mat_signature_binary, 
              file.path(dir_res_viz, cfm, 'mat_cell_is_using_mm.rds'))

    #------ Heatmap mat_patient_cellpct with patient groups -----
    mat_viz <- t(mat_patient_cellpct)
    col_pal <- colorRamp2(
        seq(from=quantile(mat_viz, .01), to=quantile(mat_viz, .99), 
            length.out = 10),
        colors = diverge_hcl(n=10, palette = 'Tropic')
        # colors = sequential_hcl(n=5, palette = 'Terrain 2')
    )
    col_anno <- HeatmapAnnotation(
        df= data.frame(
            PCR_status=patient_info_pcr[colnames(mat_viz)],
            NMF4=patient_info_nmf4[colnames(mat_viz)]), 
        col = list(PCR_status=pal_pcr, NMF4=pal_nmf4), which = 'column')
    
    for (ii in c(0, 1, 2, 3, 4)) {
        ## 0: clustering and coloring are independent
        ## 1: coloring prioritize pCR/RD info then NMF4. Clustering is done per slice.
        ## 2: coloring prioritize NMF4 then pCR/RD.  Clustering is done per slice.
        print(ii)
        col_split <- NULL
        if (ii == 0) {col_split <- NULL} 
        if (ii == 1) {col_split <- paste0(
            patient_info_pcr[colnames(mat_viz)],
            patient_info_nmf4[colnames(mat_viz)], sep='_')}
        if (ii == 2) {col_split <- paste0(
            patient_info_nmf4[colnames(mat_viz)], 
            patient_info_pcr[colnames(mat_viz)], sep='_')}
        if (ii == 3) {col_split <- patient_info_pcr[colnames(mat_viz)]}
        if (ii == 4) {col_split <- patient_info_nmf4[colnames(mat_viz)]}
        
        p <- Heatmap(
            mat_viz, name = 'fraction',
            col = col_pal, 
            show_row_names = TRUE,
            show_column_names = TRUE,
            top_annotation = col_anno,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 4),
            column_split = col_split,
            cluster_columns = ifelse(ii == 0, TRUE, FALSE), 
            clustering_method_columns = 'ward.D2',
            cluster_column_slices = T,
            column_title = NULL,
            border = T,
            row_dend_width  = unit(5, 'mm'),
            column_dend_height  = unit(5, 'mm'),
            row_gap = unit(0, 'cm'), column_gap = unit(0, 'cm'),
            # use_raster=T, raster_device = 'CairoPNG'
        )
        pdf(file.path(dir_res_viz, cfm, 
                      sprintf('heatmap.cell_fraction.%s.pdf', ii)), 
            width = 7.5, height = 4, useDingbats = F)
        draw(p)
        dev.off()
    }    

    #------ Violin plot to find which MM showed diff fraction ------
    # - pCR/RD
    # - NMF4
    patient_info <- NULL
    mm_avail <- setdiff(sort(unique(df_patient_cellpct$signature)), 'Unresolved'); print(mm_avail)
    for (z in c('nmf4', 'pcr')) {
        cat(z, '...\n')
        if (z == 'pcr') { patient_info <- patient_info_pcr; pal_use <- pal_pcr}
        if (z == 'nmf4') {patient_info <- patient_info_nmf4; pal_use <- pal_nmf4}
        print(table(patient_info))
        df_patient_cellpct$grouping <- patient_info[df_patient_cellpct$patient]
        # View(df_patient_cellpct)
        
        # Violin plot for each meta-module
        signature_names <- colnames(mat_patient_cellpct)
        p_list <- list()
        for (x in signature_names) {
            # next()
            # cat(x, ' ')
            tmp <- df_patient_cellpct %>% filter(grouping!='Unknown', signature==x)
            tmp_xlab <- tmp %>% dplyr::group_by(grouping) %>% 
                dplyr::summarise_at('pct', mean) %>% deframe() %>%
                formatC(digits = 2) %>% pretty_table2str()
            if (z=='pcr'){tmp_xlab <- tmp_xlab[c('pCR', 'RD')]}
            p <-  tmp %>%
              ggplot(aes(x=grouping, y=pct)) +
              geom_boxplot(outlier.shape=NA) +
              geom_quasirandom(aes(color=grouping), size=2, alpha=0.5, pch=16) +
              
              # geom_violin(alpha=0, width=.5) +
              ggpubr::stat_mean(pch=5, cex=2, color='gold') 
            
            p <- p + labs(x='', y='%cells using the module', 
                          title=str_remove_all(x, 'HALLMARK_') %>% str_remove('Signature_') %>% str_remove('_UCell')) +
              scale_x_discrete(limits=names(tmp_xlab), labels=tmp_xlab) +
              theme_pubr(base_size = 6) + 
              scale_y_continuous(labels=percent) + 
              coord_cartesian(ylim=c(0,1)) +
              scale_color_manual(values=pal_use) +
              rremove('legend') +
              theme(axis.ticks.length=unit(0.07,"inch"))
            if (z == 'pcr') {p <- p+stat_compare_means(method = "wilcox.test", aes(label = paste0("p = ", after_stat(p.format)))); p_width <- 2}
            if (z == 'nmf4'){p <- p+stat_compare_means(aes(label = paste0("p = ", after_stat(p.format)))); p_width <- 4}
            # p
            rm(tmp_xlab)
            ggsave(file.path(dir_res_viz, cfm, 
                             sprintf('VlnPlot.by%s.%s.cell_pct.pdf', z, x) ),
                   p, width = p_width, height = 3)
            p_list<- c(p_list, list(p))
        }; cat('.\n')
        ggsave(file.path(dir_res_viz, cfm, 
                         sprintf('VlnPlot.by%s.%s.cell_pct.pdf', z, 'compact') ),
               ggarrange(plotlist = p_list, nrow=1), width = p_width * length(p_list), height = 3, limitsize = F)
        
        # Heatmap summarizing fraction
        
        p <- df_patient_cellpct %>% 
            dplyr::group_by_at(c('grouping', 'signature')) %>% 
            dplyr::summarise(pct = mean(pct)) %>%
            ggplot(aes_string(x='grouping', y='signature')) +
            geom_tile(aes_string(fill='pct'), color='black')  +
            scale_y_discrete(limits=rev(c(mm_avail, 'Unresolved')))  +
            geom_text(aes(label=scales::percent(pct, accuracy = 0.01))) +
            scale_fill_binned_sequential() +
            theme_pubr(legend = 'right') + rremove('x.axis') + rremove('y.axis')
        print(p)
        if (z == 'pcr') {
            p <- p + scale_x_discrete(limits = c('pCR', 'RD'))
        }
        ggsave(file.path(dir_res_viz, cfm,
                         sprintf('heatmap.cell_fraction.%s.pdf', z)),
               p, width = ifelse(z=='pcr', 4, 6), height = 6,useDingbats = F)
        
        
        
        cat('Done', z, '\n')
    }
    
    #------ Statistically test meta-module pCR diff using all patients ------
        
    p <- df_patient_cellpct %>% 
        dplyr::filter(PCR %in% c('pCR', 'RD')) %>%
        dplyr::group_by_at(c('PCR', 'signature')) %>% 
        dplyr::summarise(pct = mean(pct)) %>%
        dplyr::group_by_at('signature') %>%
        dplyr::mutate(diff_fraction = pct[PCR=='RD'] - pct[PCR=='pCR']) %>%
        dplyr::select(signature, diff_fraction) %>% unique() %>%
        dplyr::mutate(color=ifelse(diff_fraction>=0, 'RD', 'pCR')) %>%
        ggplot(aes(y=signature, x=diff_fraction)) +
        geom_col(aes(fill=color)) +
        scale_y_discrete(limits=rev(c(mm_avail, 'Unresolved')))  +
        scale_x_continuous(labels=scales::percent_format(accuracy = 1)) +
        geom_vline(xintercept = 0, lty='dashed') +
        scale_fill_manual(values=pal_pcr) + rremove('legend')
    ggsave(file.path(dir_res_viz, cfm, sprintf('lolipop.cell_fraction_diff.PCR.pdf')),
           p, width = 4, height = 6,useDingbats = F)
    
    p <- df_patient_cellpct %>% 
        dplyr::filter(PCR %in% c('pCR', 'RD')) %>%
        dplyr::group_by_at(c('PCR', 'signature')) %>% 
        dplyr::summarise(pct = mean(pct)) %>%
        dplyr::group_by_at('signature') %>%
        dplyr::mutate(diff_fraction = pct[PCR=='RD'] - pct[PCR=='pCR']) %>%
        dplyr::select(signature, diff_fraction) %>% unique() %>%
        dplyr::mutate(color=ifelse(diff_fraction>=0, 'RD', 'pCR')) %>%
        dplyr::mutate(abs_diff_fraction = abs(diff_fraction)) %>%
        ggplot(aes(y=signature, x=abs_diff_fraction)) +
        geom_col(aes(fill=color)) +
        scale_y_discrete(limits=rev(c(mm_avail, 'Unresolved')))  +
        scale_x_continuous(labels=scales::percent_format(accuracy = 1)) +
        scale_fill_manual(values=pal_pcr) + rremove('legend') +
        theme(legend.position = 'right')
    ggsave(file.path(dir_res_viz, cfm, sprintf('lolipop.cell_fraction_diff.PCR.colored.pdf')),
           p, width = 5, height = 6,useDingbats = F)
    
    
    mm_viz <- unique(df_patient_cellpct$signature)
    
    #------ Statistically test meta-module pCR diff using each NMF group ------
    for (z in c(paste0('ARC', 1:4))) {
        
        p <- df_patient_cellpct %>% 
            dplyr::filter(NMF4 == z) %>%
            dplyr::filter(PCR %in% c('pCR', 'RD')) %>%
            dplyr::group_by_at(c('PCR', 'signature')) %>% 
            dplyr::summarise(pct = mean(pct)) %>%
            dplyr::group_by_at('signature') %>%
            dplyr::mutate(diff_fraction = pct[PCR=='RD'] - pct[PCR=='pCR']) %>%
            dplyr::select(signature, diff_fraction) %>% unique() %>%
            dplyr::mutate(color=ifelse(diff_fraction>=0, 'RD', 'pCR')) %>%
            ggplot(aes(y=signature, x=diff_fraction)) +
            geom_col(aes(fill=color)) +
            scale_y_discrete(limits=rev(c(mm_avail, 'Unresolved')))  +
            scale_x_continuous(labels=scales::percent_format(accuracy = 1)) +
            geom_vline(xintercept = 0, lty='dashed') +
            scale_fill_manual(values=pal_pcr) + rremove('legend') +
            labs(subtitle = z)
        ggsave(file.path(dir_res_viz, cfm, sprintf('lolipop.cell_fraction_diff.%s_only.PCR.pdf', z)),
               p, width = 4, height = 6,useDingbats = F)
        
        mm_viz <- unique(df_patient_cellpct$signature)


        
        
    }
    cat('Done', cfm, '\n\n')    
}


#-------------------------- Fraction diff for the NMF4 --------------------------    

#--------------------------
# Strategy1: 
# step1. kruskal-wallis test to see if the four groups displayed sig diff (FDR qval)
# step2. for each signature, pairwise.wilcox.test to see which two groups diff sig
#--------------------------
for (cfm in classifier_method) {
    fs::dir_create(file.path(dir_res_viz, cfm))
    cat(cfm, '\n')
    df_patient_cellpct <- read_rds( file.path(dir_res_viz, cfm, 'df_patient_cellpct.rds') )
    mat_patient_cellpct <- read_rds( file.path(dir_res_viz, cfm, 'mat_patient_cellpct.rds') )
    # Heatmap(mat_patient_cellpct)
    colnames(mat_patient_cellpct)
    rownames(mat_patient_cellpct)
    
    mat_group_cellpct <- ruok::mat_tapply(
        mat_patient_cellpct, 
        patient_info_nmf4[rownames(mat_patient_cellpct)], FUN = mean)
    
    dim(mat_patient_cellpct)
    head(df_patient_cellpct)
    mm_viz <- setdiff(unique(df_patient_cellpct$signature), 'Unresolved')
    print(mm_viz)
    
    # step1
    mm_kruskal_pval <- sapply(seq_along(mm_viz), function(m) {
        data_s <- df_patient_cellpct[df_patient_cellpct$signature == mm_viz[m], ]
        boxplot(pct ~ NMF4, data = data_s, ylim=c(0,1), col=pal_nmf4, main=mm_viz[m])
        ks_obj <- kruskal.test(pct ~ NMF4, data = data_s)
        ks_obj$p.value
    })    
    names(mm_kruskal_pval) <- mm_viz
    mm_kruskal_pval <- mm_kruskal_pval[!is.nan(mm_kruskal_pval)]
    mm_kruskal_qval <- p.adjust(mm_kruskal_pval, method = 'fdr')
    print(table(mm_kruskal_pval < 0.05))
    print(table(mm_kruskal_qval < 0.05))
    print(cbind(mm_kruskal_pval, mm_kruskal_qval))
    
    mm_sig <- names(mm_kruskal_qval)[mm_kruskal_qval < 0.05]
    print(mm_sig)
    mm_sig <- names(mm_kruskal_qval)
    # step2
    # pairwx_obj <- pairwise.wilcox.test(data_s$pct, data_s$NMF4,
    #                                    p.adjust.method = "fdr", 
    #                                    correct=FALSE)
    # as.data.frame(pairwx_obj$p.value)
    for (m in seq_along(mm_sig)) {
        cat(mm_sig[m])
        data_s <- df_patient_cellpct[df_patient_cellpct$signature == mm_sig[m], ]
        pairwx_pval <- compare_means(
            pct ~ NMF4, data_s, 
            p.adjust.method = 'fdr', 
            paired = F, 
            method = "wilcox.test",
            symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                               symbols = c("***", "**", "*", "ns")))
        
        highest_group <- names(which.max(mat_group_cellpct[, mm_sig[m]]))
        pairwx_pval <- pairwx_pval[with(pairwx_pval, group1==highest_group | group2==highest_group), ]
        pairwx_pval$y.position <- c(0.75, 0.85, 0.95) ## fixed for visualization purpose
        p <- ggplot(data_s, aes(x=NMF4, y=pct)) +
            ggbeeswarm::geom_beeswarm(aes(color=NMF4), cex=2) +
            # geom_violin(alpha=0, width=.5) +
            geom_boxplot(alpha=0, width=.5) +
            ggpubr::stat_mean(pch=5, cex=2) +
            stat_pvalue_manual(pairwx_pval, label='p.adj') +
            coord_cartesian(ylim=c(0, 1)) + 
            scale_y_continuous(labels=percent) +
            scale_color_manual(values=pal_nmf4) +
            labs(title = mm_sig[m], y = '% cells using the signature',
                 caption = sprintf('Kruskal-Wallis test qval=%s', 
                                   signif(mm_kruskal_qval[mm_sig[m]], digits = 3))) +
            rremove('legend') + rremove('x.title')
        ggsave(filename = file.path(
            dir_res_viz, cfm, 
            sprintf('VlnPlot.manual_test.bynmf4.%s.cell_pct.%s.pdf', mm_sig[m], 
                    ifelse(mm_kruskal_qval[mm_sig[m]]<0.05, 'is_sig', 'no_sig'))),
            plot = p, width = 5, height = 3.5, useDingbats = F)
    }
    
}