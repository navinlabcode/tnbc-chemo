library(Seurat)
library(tidyverse)
library(readr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(ggsignif)
source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
theme_set(theme_pubr(base_size = 6))
## First std.AddModuleScore.R to create the 'addmodulescore.csv'
## and to attach the data frames.
## Compared the archetypes across the three lineages

#------ tumor cells ------
dir_res <- file.path(
    "/volumes/USR1/yyan/project/tnbc_pre_atlas",
    "rds_rna-integrate/pat102/lv01.aneuploidy_tri_type.aneuploid.pure5"
)
fpath_sr3_ready <- file.path(dir_res, "ready.sr3.rds")
sr3 <- read_rds(fpath_sr3_ready)
df_cellmeta <- read_rds(file.path(dir_res, "sr3_metadata.df.rds"))

## new module scores of siyuan's HBCA signature
if (T) {
    ## AUCell
    df_signature_res <- read_rds(
        file.path(dir_res, "modulescore_hbca_epithelial_lineage_siyuan_head15", "aucellscore.df.rds")
    )
    dir_res_viz <- file.path(dir_res, "modulescore_hbca_epithelial_lineage_siyuan_head15", "viz_aucell")
    fs::dir_create(dir_res_viz)
    ## Module score
    df_signature_res <- read_rds(
        file.path(dir_res, "modulescore_hbca_epithelial_lineage_siyuan_head15", "modulescore.df.rds")
    )
    dir_res_viz <- file.path(dir_res, "modulescore_hbca_epithelial_lineage_siyuan_head15", "viz_module_score")
    fs::dir_create(dir_res_viz)
}


all(rownames(df_signature_res) %in% Cells(sr3))


module_names <- colnames(df_signature_res)
print(head(module_names))
sr3 <- AddMetaData(sr3, metadata = df_signature_res)

df <- FetchData(sr3, c(module_names, "patient", "fNMF"))
head(df)
if (F) {
    ## module score
    df_viz <- df %>%
        dplyr::group_by(patient) %>%
        dplyr::summarise(
            grouping = unique(fNMF),
            nCells = n(),
            Basal = mean(module_score_basal, rm.na = T),
            LumSec = mean(module_score_lumsec, rm.na = T),
            LumHR = mean(module_score_lumhr, rm.na = T)
        )
}
if (F) {
    ## UCell
    df_viz <- df %>%
        dplyr::group_by(patient) %>%
        dplyr::summarise(
            grouping = unique(fNMF),
            nCells = n(),
            Basal = mean(Signature_Basal_UCell, rm.na = T),
            LumSec = mean(Signature_LumSec_UCell, rm.na = T),
            LumHR = mean(Signature_LumHR_UCell, rm.na = T)
        )
}
if (T) {
    ## AUCell score
    df_viz <- df %>%
        dplyr::group_by(patient) %>%
        dplyr::summarise(
            grouping = unique(fNMF),
            nCells = n(),
            Basal = mean(aucell_score_basal, rm.na = T),
            LumSec = mean(aucell_score_lumsec, rm.na = T),
            LumHR = mean(aucell_score_lumhr, rm.na = T)
        )
}
df_viz$log10nCells <- log10(df_viz$nCells)

df_viz$grouping <- factor(df_viz$grouping, levels = get_levels(sort(unique(df_viz$grouping))))
head(df_viz)

library(tidyr)

df_viz_long <- df_viz %>%
    pivot_longer(
        cols = c("Basal", "LumSec", "LumHR"),
        names_to = "lineage",
        values_to = "score"
    )
head(df_viz_long)
# view(df_viz_long)
df_viz_long$lineage <- factor(as.character(df_viz_long$lineage),
    levels = c("LumSec", "LumHR", "Basal")
)
df_viz_long %>%
    dplyr::group_by(grouping, lineage) %>%
    dplyr::summarise(score = mean(score))



pal_hbca_epi <- c("LumSec" = "royalblue", "LumHR" = "deeppink", "Basal" = "seagreen")

library(ggbeeswarm)

df_viz_long %>%
    ggplot(aes(x = grouping, y = score, fill = lineage)) +
    # geom_boxplot(outlier.shape=NA, position = position_dodge(width = .8)) +
    geom_quasirandom(
        aes(size = nCells),
        alpha = 1, pch = 21, dodge.width = .8
    ) +
    # ggpubr::stat_mean(pch=18, cex=3, color="gold",
    #                   position = position_dodge(width = .8)) +
    stat_summary(
        position = position_dodge(width = .8),
        fun.y = "mean", geom = "point"
    ) +
    stat_summary(
        position = position_dodge(width = .8),
        fun.data = "mean_se", geom = "errorbar",
        color = "black"
    ) +
    ggpubr::stat_compare_means(label = "p.format") +
    scale_fill_manual(values = pal_hbca_epi)



p <- df_viz_long %>%
    ggplot(aes(x = lineage, y = score)) +
    facet_wrap(~grouping, nrow = 1) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_quasirandom(
        # aes(size = nCells),
        size = 2,
        alpha = 1, pch = 21, fill = "lightgrey"
    ) +
    stat_summary(fun.y = "mean", geom = "point", fill = "gold", pch = 23, size = 3) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", color='black') +
    # ggpubr::stat_compare_means(label = "p.format") +
    ggpubr::stat_compare_means(
        label = "p.format",
        comparisons = list(
            c("LumSec", "LumHR"),
            c("LumSec", "Basal"),
            c("LumHR", "Basal")
        )
    ) +
    # scale_fill_manual(values = pal_hbca_epi) +
    scale_x_discrete(limits = names(pal_hbca_epi)) +
    theme(axis.ticks.length = unit(0.1, "inch"))

# p
ggsave(file.path(dir_res_viz, "boxplot.module_score.lineage_vs_NMF4.with_test.pdf"),
    p,
    width = 6, height = 3.5
)

p <- df_viz_long %>%
    ggplot(aes(x = lineage, y = score)) +
    facet_wrap(~grouping, nrow = 1) +
    # geom_boxplot(outlier.shape=NA) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_quasirandom(
        size = 2, alpha = 1, pch = 21, fill = "lightgrey"
    ) +
    stat_summary(fun.y = "mean", geom = "point", fill = "gold", pch = 23, size = 2) +
    # stat_summary(fun.data = "mean_se", geom = "errorbar", color='black') +
    ggpubr::stat_compare_means(label = "p.format") +
    # scale_fill_manual(values = pal_hbca_epi) +
    scale_x_discrete(limits = names(pal_hbca_epi)) +
    theme(axis.ticks.length = unit(0.1, "inch"))

# p
ggsave(file.path(dir_res_viz, "boxplot.module_score.lineage_vs_NMF4.with_ks_test.pdf"),
    p,
    width = 6, height = 3.5
)

head(df_viz_long)

p <- ggboxplot(
    df_viz_long,
    x = "lineage", y = "score", color = "grouping", palette = pal_nmf4,
    outlier.shape = NA, fill = NA
) + geom_quasirandom(aes(color = grouping), size = .4, dodge.width = .8)

if (F) {
    stat.test <- df_viz_long %>%
        group_by(lineage) %>%
        rstatix::wilcox_test(score ~ grouping, p.adjust.method = "fdr")
    head(stat.test)
    stat.test$p.adj.txt <- signif(stat.test$p.adj, digits = 3)
    stat.test$p.adj.txt[stat.test$p.adj >= 0.05] <- "ns"
    # stat.test$p.adj.signif[stat.test$p.adj.signif == "ns"] <- "" # remove ns
}
if (T) {
    library(rstatix)
    stat.test <- df_viz_long %>%
        group_by(lineage) %>%
        dunn_test(., score ~ grouping, p.adjust.method = "BH", detailed = FALSE)
    stat.test$p.adj.txt <- signif(stat.test$p.adj, digits = 3)
    stat.test$p.adj.txt[stat.test$p.adj >= 0.05] <- "ns"
}
stat.test <- stat.test %>%
    add_xy_position(x = "lineage", dodge = 0.8)
p <- p + stat_pvalue_manual(
    stat.test,
    label = "p.adj.txt", tip.length = 0.007,
    bracket.nudge.y = 0.5
)
ggsave(file.path(dir_res_viz, "boxplot.module_score.NMF4_vs_lineage.with_test.2.pdf"),
    p,
    width = 0.8*3+1, height = 5, useDingbats = F
)
