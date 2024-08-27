# library(tidyverse)
#--------------------------
# Any gene signature score or continuous variable ~ OS in METARBRIC
# - archetypes
# - gene-based machine learning
#
# Yun Yan
#--------------------------
library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(ruok)
library(forcats)

dir_proj <- file.path(
    '.','survival_analysis', 
    'metabric_tnbc_chemoyes')
#-------------------------- Libs --------------------------    
dir_lib <- './other/METABRIC/'
lib_expression <- file.path(dir_lib, 'tnbc_chemoYES_expression.mat.rds') 
lib_survival <- file.path(dir_lib, 'tnbc_chemoYES_survival.df.rds') 
lib_pat_names <- file.path(dir_lib, 'tnbc_chemoYES.pat_names.txt') 


#------ Prepare data ------

library(readr)
if (! file.exists(lib_survival)) {
    if (T) {
        metabric1 <- read_delim("./other/METABRIC/table_S2_revised.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
        metabric2 <- read_delim("./other/METABRIC/table_S3_revised.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
        dim(metabric1)
        dim(metabric2)
        identical( colnames(metabric1), colnames(metabric2) )
        metabric <- rbind(metabric1, metabric2)
        metabric$OS_STATUS <- ruok::replace_vector(
            metabric$last_follow_up_status,
            c('a' = 0,
              'd' = 1, 
              'd-d.s.' = 1, 
              'd-o.c.' = 2)
        ); metabric$OS_STATUS <- as.numeric(metabric$OS_STATUS)
        table(metabric$OS_STATUS, useNA='ifany')
        # 0    1    2 <NA> 
        #     1081  625  266   20
        metabric$OS_DAYS <-  metabric$T
        metabric$OS_MONTHS <-  metabric$OS_DAYS / 30
        
        data_expression <- read_rds(file.path(dir_lib, 'data_expression.df.rds'))
        data_pat_names <- colnames(data_expression)[c(-1, -2)]
    }

    colnames(metabric)
    table(metabric$ER.Expr)
    sum(grepl('CT', metabric$Treatment))
    tnbc <- metabric %>% 
        dplyr::filter(ER.Expr == '-', 
                      Her2.Expr == '-', 
                      PR.Expr == '-')
    sum(grepl('CT', tnbc$Treatment))
    table(tnbc$Treatment) 
    tnbc_chemo <- tnbc[grepl('CT', tnbc$Treatment), ]
    table(tnbc_chemo$ER_IHC_status, tnbc_chemo$ER.Expr)
    write_rds(tnbc_chemo, lib_survival)
    yun_obs_names <- tnbc_chemo$METABRIC_ID
    write_lines(yun_obs_names, lib_pat_names)
}
if (! file.exists(lib_expression)) {
    yun_obs_names <- read_lines('.other//METABRIC/tnbc_chemoYES.pat_names.txt')
    str(yun_obs_names) # 167
    
    data_expression <- read_delim(".other//METABRIC/data_expression.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
    data_expression <- read_rds(file.path(dir_lib, 'data_expression.df.rds'))
    print(dim(data_expression))
    str(intersect(yun_obs_names, colnames(data_expression)))
    mat_expression <- data_expression[, intersect(yun_obs_names, colnames(data_expression)) ]
    mat_expression <- as.matrix(mat_expression)
    rownames(mat_expression) <- data_expression$Hugo_Symbol
    dim(mat_expression)
    write_rds(data_expression, file.path(dir_lib, 'data_expression.df.rds'))
    write_rds(mat_expression, lib_expression)
}

df_surv <- read_rds(lib_survival); nrow(df_surv)
try(df_surv$Patient.ID <- df_surv$METABRIC_ID)
mat <- read_rds(lib_expression); ncol(mat)  ## log-transformed
shared_pat <- intersect( df_surv$Patient.ID, colnames(mat) ); str(shared_pat)
df_surv <- df_surv[match(shared_pat, df_surv$Patient.ID), ]
mat <- mat[, match(shared_pat, colnames(mat))]
table(df_surv$OS_STATUS)
nrow(df_surv)
table(df_surv$Relapse.Free.Status, useNA='ifany')
table(df_surv$Overall.Survival.Status, useNA='ifany')
view(df_surv)
stopifnot( identical(df_surv$Patient.ID, colnames(mat)) )
#-------------------------- Inputs --------------------------    
# HLA
if (F) {
    dir_res <- './survival_analysis/survival_metabric_chemoYES_HLA'
    dir.create(dir_res)
    setwd(dir_res)
    module_content <- list(
        HLA_Class1 = c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G'),
        HLA_Class2 = c('HLA-DPA1', 'HLA-DPB1', 
                       'HLA-DQA1', 'HLA-DQB1',
                       'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
                       'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB')
    )
    module_content_use <- module_content
}

## metaprograms
if (F) {
    module_content <- read_rds(
        file.path('.',
                  'metamodule_fnmf', 
                  'MM_alt_clean_byscore_wardD2', 'deliver.mm_markers.rds')
    )
    module_content <- read_rds(
        file.path('.',
                  'metamodule_fnmf', 
                  'MM_alt_clean_byscore_wardD2', 'deliver.mm_top_genes.rds')
    )
    str(module_content)
    module_focuse <- setdiff(names(module_content), 'M14')
    module_focuse_bio <- c('G2/M', 'Mitochondria', 'Ribosome', 'Stress', 'Interferon',
                           'HLA','S/G1', 'Hypoxia', 'Basal','EMT',
                           'LumSec', 'Cholestero', 'StressER')
    names(module_focuse_bio) <- module_focuse
    publish_top_markers <- 30
    module_content_use <- lapply(module_content, function(x) head(x, publish_top_markers))
    module_content_use <- module_content_use[module_focuse]
    str(module_content_use)
    names(module_content_use) <- paste0(
        names(module_content_use), '__',
        make.names(module_focuse_bio[names(module_content_use)]))
    db <- module_content_use
    
    dir_res <- file.path(dir_proj, 'metamodule_fnmf', 
                         'logfpkm.Gene_marker.Method_ssgsea')
    fs::dir_create(dir_res)
    wd0 <- getwd(); print(wd0)
    setwd(dir_res)    
    
}


#------ psbulk fNMF4 ------
if (TRUE) {
    nmf_genes <- read_rds('./psbulk/fastnmf/rank4/deliver.nmf_markers.top.rds')
    nmf_genes <- lapply(nmf_genes, head, n=30)
    str(nmf_genes)
    dir_res <- file.path(dir_proj, 'psbulk_nmf4', 'Gene_top30.Method_ssgsea') # use in paper
    fs::dir_create(dir_res)
    setwd(dir_res)
    module_content_use <- nmf_genes
    if (F){
        library(fgsea)
        vanderbilt_genes <- gmtPathways('./other/vanderbilt_tnbc/vanderbilt_tnbc_subtypes_gene_signatures_plus.gmt')
        str(vanderbilt_genes)
        module_content_use <- c(nmf_genes, vanderbilt_genes)
    }
    
}

fs::dir_create(dir_res)
setwd(dir_res)
str(db)
#-------------------------- GSVA --------------------------    

# view(df_surv)
library(GSVA); library(GSEABase)
library(msigdbr)

db_msigdb_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol)
db_msigdb_H <- deframe_to_list(as.data.frame(db_msigdb_H))
names(db_msigdb_H) <- make.names(names(db_msigdb_H))


gene_set_gsva <- c(db)
str(gene_set_gsva)

gsva_method_choices <- c('gsva', 'ssgsea', 'expr')
gsva_method <- 'ssgsea'

print(gsva_method)
#------ Calculate/Load the gsva matrix ------
if (gsva_method == 'gsva') {
    idx <- apply(mat, 1, function(x) var(x, na.rm = TRUE) != 0); print(table(idx))
} else {
    idx <- 1:nrow(mat)
}
if (file.exists(file.path(dir_res, sprintf('gsva_output.%s.mat.rds', gsva_method)))){
    gsva_mat <- read_rds(file.path(dir_res, sprintf('gsva_output.%s.mat.rds', gsva_method)))
} else {
    cat('calculating', gsva_method, '...')
    if (gsva_method %in% c('gsva', 'ssgsea')) {
        gsva_mat <- gsva(
            expr = mat[idx, ],
            gset.idx.list = gene_set_gsva, 
            method=gsva_method,
            min.sz = 1, max.sz=1000,
            parallel.sz=40)
    } else {
        ## gene expression only
        str(sapply(gene_set_gsva, function(gn) intersect(gn, rownames(mat))))
        gsva_mat <- sapply(gene_set_gsva, function(gn){
            log(colSums(exp(mat[intersect(gn, rownames(mat)), ])))
        })
        gsva_mat <- t(gsva_mat)
    }    
    write_rds(
        gsva_mat, 
        file.path(dir_res, sprintf('gsva_output.%s.mat.rds', gsva_method)) )
}
print(dim(gsva_mat)) ## signatures x patients
head(colnames(gsva_mat))
print(rownames(gsva_mat))

zscore_clip <- function(x, zmax=3, zmin=-3){
    pmin(pmax(x, zmin), zmax)
}
for (gs_x in rownames(gsva_mat) ) {
    pdf(file.path(dir_res, sprintf('%s_%s_score_distribution.pdf', gs_x, gsva_method)), useDingbats = F)
    hist(gsva_mat[gs_x, ], breaks = 30, main=gs_x, xlab=gsva_method)
    abline(v=c(-.1, .1), col='red')
    dev.off()
    if (gsva_method %in% c('ssgsea', 'expr')) {
        pdf(file.path(dir_res, sprintf('%s_%s_zscored_distribution.pdf', gs_x, gsva_method)), useDingbats = F, onefile = T)
        hist(zscore_clip(scale(gsva_mat[gs_x, ])), breaks = 30, main=gs_x, xlab=sprintf('%s clipped zscore', gsva_method))
        abline(v=c(-.1, .1), col='red')
        hist(scale(gsva_mat[gs_x, ]), breaks = 30, main=gs_x, xlab=sprintf('%s zscore', gsva_method))
        abline(v=c(-.1, .1), col='red')
        dev.off()
    }
}
names(gene_set_gsva)

#--------------------------
#---- Cox proportional-hazards regression ----
# - gsva: the original gsva score
# - ssgsea: zscored ssgesa
# - expr: zscored log(sum FPKM)
#--------------------------
library("survival")
library("survminer")

df <- df_surv
df <- as.data.frame(df)
#------ Merge the clinical data frame and the features ------
df$Overall.Survival.Status
df$Overall.Survival..Months.
df$Relapse.Free.Status
df$Relapse.Free.Status..Months.
colnames(df)

if (T) {
    # official metabric
    head(df$OS_DAYS)
    cat_y_choices <- c('OS_STATUS')
    cat_x_choices <- c('OS_DAYS')
    df_cox <- df[, c('Patient.ID', 'age_at_diagnosis', 
                     'NPI', 'stage', 
                     cat_y_choices, cat_x_choices)]
}
dim(gsva_mat)
gsva_mat_use <- t(gsva_mat)[as.character(df_cox$Patient.ID), , drop=F]
gsva_mat_use[1:3, 1:3]

library(forcats)
if (gsva_method %in% c('ssgsea', 'expr')) {
    cat('scale...')
    gsva_mat_zscore <- scale(gsva_mat_use)
    df_cox <- cbind(df_cox, as.data.frame(gsva_mat_zscore))
} else {
    df_cox <- cbind(df_cox, as.data.frame(gsva_mat_use))
}

gs_names <- colnames(gsva_mat_use); str(gs_names)
colnames(df_cox)

#------ muti-variate ------
    
cat_i <- 1
for (cat_i in seq_along(cat_y_choices)) {
    df_cox_i <- df_cox
     
    cat(cat_i)
    df_cox_i$time <- df_cox_i[, cat_x_choices[cat_i]]
    df_cox_i$status <- df_cox_i[, cat_y_choices[cat_i]]
    print(table(df_cox_i$status, useNA='ifany'))
    df_cox_i <- df_cox_i[df_cox_i$status %in% c(0, 1), ]
    print(table(df_cox_i$status, useNA='ifany'))

    if (T) {
        res <- coxph(
            as.formula(paste0(
                'Surv(time, status) ~ ', 
                paste(c('age_at_diagnosis', 
                        'NPI',
                        # 'stage',
                        gs_names), 
                      collapse = '+'))), 
            data =  df_cox_i)
    }
    print(summary(res))

    library(broom)
    p <- try(ggforest(res))
    p
    if ('try-error' %in% class(p)) {cat('[warn]', cat_y_choices[cat_i], ' did not work!!!'); next()}
    ggsave(file.path(dir_res, sprintf('ggforest.%s.%s.pdf', 
                                      gsva_method,
                                      cat_y_choices[cat_i])),
           plot = p,
           width = 9, height = 7)
}

