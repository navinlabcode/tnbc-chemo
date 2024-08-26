#--------------------------
# plot for the cell cycle phase comparison of pcr v.s. rd
#
# Yun Yan
#--------------------------
dir_res <- '.'
fpath_sr3_ready <- file.path(dir_res, 'ready.sr3.rds')
# sr3 <- read_rds(fpath_sr3_ready)
df_cellmeta <- read_rds(file.path(dir_res, 'sr3_metadata.df.rds'))

dir_res_viz <- './viz_signature_sr3'
fs::dir_create(dir_res_viz)

#-------------------------- count cell freq
colnames(df_cellmeta)
df_p1 <- df_cellmeta %>% dplyr::group_by(patient) %>%
    dplyr::summarise(nCell = n(), 
                     nS = sum(Phase == 'S'),
                     nG1 = sum(Phase == 'G1'),
                     nG2M = sum(Phase == 'G2M'),
                     PCR_status = unique(PCR_status)) %>%
    dplyr::mutate(pS = nS / nCell,
                  pG1 = nG1 / nCell,
                  pG2M = nG2M / nCell)

head(df_p1)
df_p1 <- df_p1 %>% dplyr::group_by(PCR_status) %>%
    dplyr::summarise(pS = mean(pS), 
                     pG1 = mean(pG1), 
                     pG2M = mean(pG2M)) %>%
    dplyr::filter(PCR_status %in% c('pCR', 'RD'))

df_p1 <- df_p1 %>% gather('color', 'value', -PCR_status)
theme_set(theme_pubr())
p1 <- ggplot(df_p1, aes(y=PCR_status, x=value, 
                        fill = color)) + 
    geom_col() +
    geom_text(aes(label=scales::percent(value)), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values=c('pG1'='#c5d642', 
                               'pS'='#f2e146', 
                               'pG2M'='#88b7c3')) +
    scale_y_discrete(limits=c('RD', 'pCR')) +
    rremove('x.axis') + rremove('x.ticks') + rremove('x.text') + rremove('y.axis')

ggsave(file.path(dir_res_viz, 'barplot.cell_cycle_proportion.PCR.pdf'),p1, width = 4, height = 3 )

#-------------------------- Cell phase score
table(df_cellmeta$Phase)
df_cellmeta

df_p2 <- df_cellmeta %>% group_by(patient, PCR_status) %>%
    summarise(S.Score = mean(S.Score), 
              G2M.Score = mean(G2M.Score)) 
pal_pcr <- c(
    'pCR' = 'seagreen3',
    'RD'='tomato1', 
    'Unknown'='black', 'Excluded'='lightgrey')
library(ggbeeswarm)
p2_a <- ggplot(df_p2, aes(x=PCR_status, y=S.Score)) +
    scale_x_discrete(limits = c('pCR', 'RD')) +
    scale_color_manual(values = pal_pcr) +
    geom_boxplot(outlier.colour = NULL, width=.5, alpha = 0) +
    geom_beeswarm(cex=2, aes(color=PCR_status)) +
    stat_compare_means() 
p2_b <- ggplot(df_p2, aes(x=PCR_status, y=G2M.Score)) +
    scale_x_discrete(limits = c('pCR', 'RD')) +
    scale_color_manual(values = pal_pcr) +
    geom_boxplot(outlier.colour = NULL, width=.5, alpha = 0) +
    geom_beeswarm(cex=2, aes(color=PCR_status)) +
    stat_compare_means() 

ggsave(file.path(dir_res_viz, 'barplot.S.Score.PCR.pdf'),p2_a, width = 3, height = 4 )
ggsave(file.path(dir_res_viz, 'barplot.G2M.Score.PCR.pdf'),p2_b, width = 3, height = 4 )

