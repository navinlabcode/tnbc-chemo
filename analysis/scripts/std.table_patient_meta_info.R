#--------------------------
# Table plot any two categorical variable of patients
# Yun Yan
#--------------------------
library(ggpubr)
library(tidyverse)
library(gtools)
library(colorspace)
source('util.ruok.R')

dir_res <-'./patient_meta_info_crosstable'
fs::dir_create(dir_res)

#------ pat artemis_id ~ subject_id ------
df_cellmeta <- read_rds('./sr3_metadata.df.rds')
pat_id_dict <- df_cellmeta[, c('patient', 'subject_id')] %>% unique() %>% tibble()
head(pat_id_dict)
#------ pcr info ------
pcr <- read_rds('~/project/tnbc_pre_atlas/meta/pCR_info_all.clean.rds') %>% tibble()
pcr$subject_id <- as.character(pcr$subject_id)
head(pcr)
sum(table(pcr$pCR_status))
table(pcr$pCR_status)
pat_pcr <- left_join(pat_id_dict, pcr, by='subject_id')
head(pat_pcr)

#--------------------------
# psbulk-NMF x pcr/non-pcr
#--------------------------

#------ nmf4 ------
psbulk_nmf <- read_rds('./psbulk/fastnmf/rank4/deliver.patient_best_nmf.rds')
sum( table(psbulk_nmf) )
head(psbulk_nmf)
cat_psbulk <- 'psbulk_nmf4'
psbulk_nmf <- enframe(psbulk_nmf, 'patient', cat_psbulk)
cat_psbulk_lv <- mixedsort(unique(unique(psbulk_nmf$psbulk_nmf4)))


pat_avail <- intersect(psbulk_nmf$patient, pat_pcr$patient); str(pat_avail)

df <- left_join(psbulk_nmf, pat_pcr, by='patient')
table(df$pCR_status)


#---------------------------
# psbulk-NMF4 x bulk vanderbilt types              ----    
#---------------------------
library(readxl)
psbulk_nmf <- read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/export_to_paper/patient_meta_info/export.sample_info.rds')
psbulk_nmf <- psbulk_nmf[, c('subject_id', 'archetype')]
cat_psbulk_lv <- sort(unique(psbulk_nmf$archetype))
bulk_vanderbilt <- readxl::read_excel('/volumes/USR1/yyan/project/tnbc_pre_atlas/meta/bulk_vanderbilt.xlsx')
# bulk_vanderbilt <- deframe(bulk_vanderbilt)
df <- inner_join(x=psbulk_nmf, y=bulk_vanderbilt)
cat_psbulk = 'archetype'

#---------------------------
# Viz Any things 
#---------------------------
df <- read_rds('/volumes/USR1/yyan/project/tnbc_pre_atlas/export_to_paper/patient_meta_info/export.sample_info.rds')

#------ viz per category ------
colnames(df)
library(patchwork)

z_opts = c('pCR_status', 'archetype')
p <- ggpubr::ggarrange(
  plotlist = lapply(z_opts, function(z) {
    qbarplot_table_cat(df[[z]], do.prop.table = F, name_x = z)
  }), 
  align = 'h', nrow=1
)
ggsave(filename = file.path(dir_res, 'barplot_overview_all_cats.pdf'), plot = p,
       width = length(z_opts) * 3 + 0.5, height = 5)

#------ visualization ------
theme_set(theme_minimal())
if (T) {
  # psbulk-nmf x pCR/RD
  a = cat_psbulk
  b = 'pCR_status'
}


df$a <- df[[a]]
df$b <- df[[b]]
mat0 <- as.matrix(table(df$a, df$b, useNA='ifany'))
rownames(mat0) <- make.names(rownames(mat0))
colnames(mat0) <- make.names(colnames(mat0))
pdf(file.path(dir_res, sprintf('%s_X_%s.raw_original.pdf', a, b)),  width = 6, height = 5, useDingbats = F)
p <- Heatmap(mat0, name = sprintf('num (N=%s)', sum(mat0)),  
             row_names_side = 'left', column_names_side = 'top',
             cluster_rows = F, cluster_columns = F, 
             col = heatmap_color_fun_cont(from=0, 
                                          to=max(mat0), 
                                          palette = 'Blues 2', rev = T),
             width = unit(4*ncol(mat0)/nrow(mat0), 'inch'),  
             height = unit(4, 'inch'), 
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (mat0[i, j] >= 0) {
                 grid.text(sprintf("%d", mat0[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
               } 
             },
             use_raster=T, raster_by_magick = TRUE)
draw(p)
dev.off()


cat_a_lvs <- sort(unique(df$a)); print(cat_a_lvs)
cat_b_lvs <- sort(unique(df$b)); print(cat_b_lvs)

p0 <- df %>% dplyr::group_by(b, a) %>% dplyr::summarise(n=n()) %>% 
    ggplot(aes(x=b, y=a, fill=n)) +
    geom_tile() +
    geom_text(aes(label=n)) +
    colorspace::scale_fill_binned_sequential() +
    labs(x=b, y=a, fill='nPatients') +
  scale_y_discrete(limits=rev(cat_a_lvs)) + 
  scale_x_discrete(limits = rev(cat_b_lvs))
p0
ggsave(file.path(dir_res, sprintf('%s_X_%s.raw.pdf', a, b)), p0, width = 5, height = 4)


#------ Focus on pCR and RD only ------
  
df2 <- df %>% 
    dplyr::filter(b %in% c('pCR', 'RD'))
p1 <- df2 %>%
    dplyr::group_by(b, a) %>% dplyr::summarise(n=n()) %>% 
    ggplot(aes(x=b, y=a, fill=n)) +
    geom_tile() +
    geom_text(aes(label=n)) +
    colorspace::scale_fill_binned_sequential() +
    labs(x=b, y=a, fill='nPatients') +
    scale_y_discrete(limits=rev(cat_a_lvs), labels=pretty_table2str(table(df2$a))) +
    scale_x_discrete(labels=pretty_table2str(table(df2$b)))
p1
ggsave(file.path(dir_res, sprintf('%s_X_%s.viz.pdf', a, b)), p1, width = 5, height = 4)


mat <- table(df2$a, df2$b)
mat <- mat[, c('pCR', 'RD')]
mat
source("util.nmf.viz.R")

pdf(file.path(dir_res, sprintf('%s_X_%s.raw.2.pdf', a, b)),  width = 6, height = 5, useDingbats = F)
p <- Heatmap(mat, name = 'N', 
        row_names_side = 'left', column_names_side = 'top',
        cluster_rows = F, cluster_columns = F, 
        col = heatmap_color_fun_cont(from=0, 
                                     to=max(mat), 
                                     palette = 'Blues 2', rev = T),
        width = unit(4*ncol(mat)/nrow(mat), 'inch'),  
        height = unit(4, 'inch'), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (mat[i, j] > 0) {
            grid.text(sprintf("%d", mat[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
          } 
        },
        use_raster=T, raster_by_magick = TRUE)
draw(p)
dev.off()

library(ggpubr)
res <- chisq.test(mat)
res$p.value
res$statistic

p2 <- res$residuals %>% as.data.frame() %>%
    ggplot(aes(y=Var1, x=Var2, color=Freq)) + 
    geom_point(size=8) + 
    colorspace::scale_color_continuous_diverging() + 
    scale_y_discrete(limits=rev(cat_a_lvs)) +
    theme_pubr(legend = 'right') +
    labs(color='pearson residuals', x=b, y=a)
p2


p3 <- res$residuals %>% as.data.frame() %>%
  ggplot(aes(x=Freq, y=Var1)) +
  geom_vline(xintercept = chi_residual_cutoff(res), lty='dashed') +
  geom_vline(xintercept = 0, lty='solid') + 
  geom_col(aes(fill=Freq>=0)) + facet_wrap(~Var2) +
  labs(x='pearson residuals', y=a, 
       caption = sprintf('chi-square (%s, %s)=%s pval=%s', 
                         res$parameter, sum(res$observed), round(res$statistic, 2), scientific(res$p.value, digits = 2))
  ) +
  scale_fill_manual(values=c(`TRUE`='orange', `FALSE`='grey'))+
  scale_y_discrete(limits=rev(cat_a_lvs)) + 
  ggthemes::theme_base() + rremove('legend') 
p3
ggsave(file.path(dir_res, sprintf('%s_X_%s.viz.chisquare_test.pdf', a, b)), p3, width = 5, height = 4)


