#--------------------------
# Form meta-modules: justify collapsing clusters of factors
# 'automate' the manual merging steps
# 
# Yun Yan (yun.yan@uth.tmc.edu)
#--------------------------

library(clustree)
library(tidyverse)
library(readr)
library(randomcoloR)
source('util.nmf.viz.R')
library(ggpubr)
library(forcats)
theme_set(theme_pubr())

dir_proj <- './metamodule_fnmf/all_ranks/'
z <- 'jaccard'
mat <- read_rds( file.path(dir_proj, 'jaccard_mat.all_ranks.in_use.rds') )
print(dim(mat))


#------ Input ------
if (T) { 
    dir_res <- file.path(dir_proj, z, 'manual_round1'); fs::dir_create(dir_res)
    nmf_mem <- read_csv(file.path(dir_proj, z, 'hclust', 
                                  'metamodule_membership.jaccard.24.hclust_ward.D2.csv'))
    nmf_mem_list <- Sys.glob( file.path(dir_proj, z, 'hclust', 'metamodule_membership*ward.D2.csv') )    
}


#------ clustree the different NMF models ------
# nmf_mem <- read_csv(nmf_mem_list[1])
head(nmf_mem$meta_module)
nmf_mem_names <- nmf_mem$index_name
str(nmf_mem_names)
nmf_mem_list <- lapply(nmf_mem_list, function(x){
    new_name <- paste0('Rank', str_extract(basename(x), pattern = '[0-9]+'))
    # metamodule_membership.jaccard.10.hclust_ward.D2.csv => 10
    y <- read_csv(x)
    y <- y[, c('index_name', 'meta_module')]
    y <- as.data.frame(y)
    rownames(y) <- y$index_name
    y <- y[nmf_mem_names, ]
    y$index_name <- NULL
    colnames(y) <- c(new_name)
    return(y)
})
library(dplyr)
df_nmf_mem <- do.call('cbind', nmf_mem_list)

p <- clustree(df_nmf_mem, prefix = 'Rank')
ggsave(file.path(dir_res, 'clustree.pdf'), p, width = 12, height = 7)

#------ inter- and intra-block similarity ------
nmf_mem <- nmf_mem[, c('index_name', 'meta_module')] %>% deframe()
table(nmf_mem)
mat <- mat[names(nmf_mem), names(nmf_mem), drop=F]
print(dim(mat))
head(nmf_mem)

# stopifnot( isSymmetric(mat) )

mat_tapply_row_col_both <- function(x, y){
    lv_opt <- gtools::mixedsort( unique(y) )
    res <- matrix(0, nrow = length(lv_opt), ncol = length(lv_opt))
    for (res_i in seq_along(lv_opt)) {
        for (res_j in seq_along(lv_opt)) {
            lv_i <- lv_opt[res_i]
            lv_j <- lv_opt[res_j]
            mat_i <- which(y == lv_i)
            mat_j <- which(y == lv_j)
            mat_v <- x[mat_i, mat_j]
            res_v <- NA
            if (lv_i == lv_j) {
                res_v <- mean(mat_v[lower.tri(mat_v)], na.rm = T)
            } else {
                res_v <- mean(mat_v, na.rm = T)
            }
            res[res_i, res_j] <- res_v
        }
    }
    rownames(res) <- colnames(res) <- lv_opt
    return(res)
}

mat_mean <- mat_tapply_row_col_both(x=mat, y=nmf_mem)

print(diag(mat_mean))
intra_block_similarity <- diag(mat_mean)


p <- enframe(intra_block_similarity) %>%
    ggplot(aes(x=fct_reorder(name, value), y=value, fill = value > 0.1)) + geom_col() +
    geom_hline(yintercept = 0.1, lty='dashed', color='grey') +
    labs(x='MM candidates', y=sprintf('intra-MM similarity (%s)', z)) +
    rremove('legend') + scale_fill_manual(values=c(`TRUE`='black', `FALSE`='red')) +
    rotate_x_text()
p
ggsave(file.path(dir_res, 'lolipop.intra_MM_candidate_similarity.pdf'), p, width = 7, height = 4)

p1 <- Heatmap(mat_mean, name = sprintf('avg %s', z),  
              cluster_rows = F, cluster_columns = F, 
              height = unit(5, 'inch'), width = unit(5, 'inch'))
p2 <- Heatmap((mat_mean > 0.1) * 1, name = sprintf('avg %s > 0.1', z), cluster_rows = F, cluster_columns = F, 
              height = unit(5, 'inch'), width = unit(5, 'inch'),
              col = c(`0` = 'grey', `1`='black'))
pdf(file.path(dir_res, sprintf('heatmap.meta_module_candidate_center_similarity.%s.pdf', z)), width = 7, height = 7, useDingbats = F)
draw(p1)
dev.off()

pdf(file.path(dir_res, sprintf('heatmap.meta_module_candidate_center_similarity_binarized.%s.pdf', z)), width = 7, height = 7, useDingbats = F)
draw(p2)
dev.off()


#------ Decisions to merge the over-clustered meta-modules ------

if (T) {
    ## For example: 
    df_nmf_mem$manual <- df_nmf_mem$Rank24
    ## Merge M01-M04 as M01
    idx <- df_nmf_mem$Rank24 %in% c('M01', 'M02', 'M03', 'M04')
    df_nmf_mem$manual[idx] <- 'M01'
}


#------ Continue ------

c(table(df_nmf_mem$manual))
library(gtools)
df_nmf_mem$manual <- factor(df_nmf_mem$manual, levels = mixedsort( unique(df_nmf_mem$manual) ))
df_nmf_mem$manual <- as.numeric(df_nmf_mem$manual)
library(stringr)
df_nmf_mem$manual <- paste0("M", str_pad(df_nmf_mem$manual, width = 2, pad = '0'))

library(tidyverse)
nmf_mem <- structure(df_nmf_mem$manual, names=rownames(df_nmf_mem))
print(enframe(c(table(nmf_mem))))
write_rds(nmf_mem, file.path(dir_res, 
                             sprintf('metamodule_membership.%s.manual.rds', z)))
stats_nmf_mem <- enframe(nmf_mem, name='index_name', value = 'meta_module') %>%
    tidyr::separate(index_name, c('sample', 'Rank', 'NMF'), sep='_', remove=F)
df <- stats_nmf_mem %>% dplyr::group_by(meta_module) %>%
    dplyr::summarise(n=length(unique(sample)))
p <- ggplot(df, aes(x=meta_module, y=n)) + 
    geom_bar(stat = 'identity') +
    geom_text(aes(label=n)) + ggpubr::rotate_x_text(90) +
    labs(y='number of sample', 
         caption = sprintf('avg %.0f samples per meta module', mean(df$n)))
p
ggsave(file.path(dir_res, sprintf('%s-%s.barplot.num_samples_per_module.pdf', 
                                  z, 'manual')), 
       p, width = 7, height = 4)
write_csv(df, file.path(dir_res, sprintf('%s-%s.num_samples_per_module.csv', 
                                         z, 'manual')))
write_csv(stats_nmf_mem, 
          file.path(dir_res, 
                    sprintf('%s-%s.metamodule_membership.csv', 
                            z, 'manual')))

#------ Heatmap ------
source('util.nmf.viz.R')
library(randomcoloR)
nmf_mem <- read_rds(file.path(dir_res, 
                              sprintf('metamodule_membership.%s.manual.rds', z)))
mat <- read_rds( file.path(dir_proj, 'jaccard_mat.all_ranks.in_use.rds') )
nmf_mem_lv <- sort(unique(nmf_mem))
nmf_mem_color <- structure(distinctColorPalette(k=length(nmf_mem_lv)), names=nmf_mem_lv)
nmf_mem_color <- structure(
    ArchR::paletteDiscrete(nmf_mem_lv),
    names=nmf_mem_lv) ; scales::show_col(nmf_mem_color)

p2 <- heatmap_cor_mat(
    cor_mat = mat[names(nmf_mem), names(nmf_mem)],
    hclust_obj = NULL,
    col = colorRamp2(
        seq(from=0, to=1, length.out = 5),
        colors = sequential_hcl(n=5, palette = 'Oslo', rev = F)
    ),
    cor_name = z, run_name = '',
    column_gap = unit(0, 'mm'), row_gap = unit(0, 'mm'), border=T,
    row_title_rot = 0, column_title_rot=90,
    row_split=nmf_mem,
    column_split=nmf_mem,
    left_annotation = rowAnnotation(MM=nmf_mem, col=list(MM=nmf_mem_color)),
    top_annotation = columnAnnotation(MM=nmf_mem, col=list(MM=nmf_mem_color)),
    topdf = file.path(dir_res,
                      sprintf('heatmap.NMFs_%s_to_modules.M%s.by_hclust_%s.pdf',
                              z, length(unique(nmf_mem)), 'manual')),
    use_raster = T, raster_device = 'CairoPNG', 
    raster_by_magick = T)

p2 <- heatmap_cor_mat(
    cor_mat = mat[names(nmf_mem), names(nmf_mem)],
    hclust_obj = NULL,
    col = colorRamp2(
        seq(from=0, to=1, length.out = 5),
        colors = sequential_hcl(n=5, palette = 'Oslo', rev = F)
    ),
    cor_name = z, run_name = '',
    column_gap = unit(0, 'mm'), row_gap = unit(0, 'mm'), border=T,
    row_title_rot = 0, column_title_rot=90,
    
    left_annotation = rowAnnotation(MM=nmf_mem, col=list(MM=nmf_mem_color)),
    top_annotation = columnAnnotation(MM=nmf_mem, col=list(MM=nmf_mem_color)),

    use_raster = T, raster_device = 'CairoPNG',
    raster_by_magick = T,
    topdf = file.path(dir_res,
                      sprintf('heatmap.NMFs_%s_to_modules.M%s.by_%s.pdf',
                              z, length(unique(nmf_mem)), 'nocut'))
)


