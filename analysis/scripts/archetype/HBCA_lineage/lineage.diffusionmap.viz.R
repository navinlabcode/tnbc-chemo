suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Project the TNBC cells onto the normal epithelial cells
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2025-05-21
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
suppressPackageStartupMessages({
    library(readr)
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    theme_set(theme_pubr(base_size = 6, legend = "right") %+replace% theme(axis.ticks.length = unit(0.1, "inch")))
    library(cli)
    library(tictoc)
    library(glue)
    library(scales)
    library(tools)
    library(fs)
    library(nanoparquet)
    library(Seurat)
    library(destiny)
    library(scattermore)
    library(ggrastr)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
})
cmdargs <- commandArgs(trailingOnly = TRUE)
if (length(cmdargs) > 0) {} else {}

cat("[done]")
timestamp()

#------ Input ------
if (F) {
    f_sr3 <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Tumor/downto_1000_by_cell_state_paper/ready.sr3.rds"
    f_dfmeta <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Tumor/downto_1000_by_cell_state_paper/sr3.metadata.df.rds"
}
if (T) {
    f_sr3 <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Tumor/downto_100_by_patient/ready.sr3.rds"
    f_dfmeta <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Tumor/downto_100_by_patient/sr3.metadata.df.rds"
}


if (F) {
    ## what if using T cells for example?
    f_sr3 <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/T/downto_1000_by_cell_state_paper/ready.sr3.rds"
    f_dfmeta <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/T/downto_1000_by_cell_state_paper/sr3.metadata.df.rds"
}

dir_res <- file.path(dirname(f_sr3), "diffusion_map_to_HBCA")
fs::dir_create(dir_res)

#------------
# The diffusion map of normal epithelial cells is kindly provided by Siyuan He.
#------------
f_ref <- "/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/cell_origin/diffusion_map/normal_epi_sub_count_dm.rds"
f_cell2celltype <- "/volumes/USR1/siyuan/DCIS_siyuan/scRNA-seq/final_round/tumor/cell_origin/diffusion_map/normal_epi_sub_anno.txt"

ref <- read_rds(f_ref)
print(class(ref))
ref_cell2celltype <- read_delim(f_cell2celltype, delim = " ", skip = 1, col_names = F)
head(ref_cell2celltype)
colnames(ref_cell2celltype) <- c("cell", "cell_type")
ref_cell2celltype <- deframe(ref_cell2celltype)
table(duplicated(names(ref_cell2celltype)))
head(ref_cell2celltype)
table(ref_cell2celltype)

df_meta <- read_rds(f_dfmeta)
colnames(df_meta)
table(df_meta$pCR_status)
table(df_meta$archetype)

sr3 <- read_rds(f_sr3)
print(sr3)


#------------------ ~~~ Predict ~~~ --------------------
library(ggnewscale)
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")

cli_h1("Predict")
query_mat <- sr3@assays$RNA@data
dim(query_mat)
query_mat <- as.matrix(query_mat)

ref_genes <- colnames(dataset(ref))
stopifnot(all(ref_genes %in% rownames(query_mat)))
print(range(dataset(ref)))
query_mat <- query_mat[ref_genes, ]
print(range(query_mat))

f_dm_predict <- file.path(dir_res, "dm_predict.rds")

if (file.exists(f_dm_predict)) {
    cli_alert_info("Loading existing prediction results from {f_dm_predict}")
    pred <- read_rds(f_dm_predict)
} else {
    cli_alert_danger("Should predict diffusion map for query matrix with {ncol(query_mat)} cells")
    # pred <- dm_predict(ref, t(query_mat))
    # write_rds(pred, f_dm_predict)
}

df_ref <- as.data.frame(ref)
df_ref <- df_ref[, setdiff(colnames(df_ref), ref_genes)]
length(colnames(df_ref))
head(rownames(df_ref)) ## cell names

ref_pstime <- rank(eigenvectors(ref)[, 1])
head(ref_pstime)
# ref_dpt <- DPT(ref)


identical(names(ref_pstime), rownames(df_ref))
df_ref$pseudotime_diffusionmap <- ref_pstime

df_pred <- as.data.frame(pred)
head(df_pred)
rownames(df_pred) <- colnames(query_mat)
identical(colnames(query_mat), Cells(sr3))

# pred_pstime <- rank(eigenvectors(pred)[, 1]); head(pred_pstime)

wanted_z_opts <- c("archetype", "pCR_status", "patient")
identical(rownames(df_pred), rownames(df_meta))
for (z in wanted_z_opts) {
    df_pred[, z] <- df_meta[, z]
}

black_cells <- !df_pred$pCR_status %in% c("pCR", "RD")
df_pred <- df_pred[!black_cells, ]

df_ref$cell_type <- ref_cell2celltype[rownames(df_ref)]

ptsize.ggscatermore <- 2
ptsize.ggrastr <- 0.2
p0 <- ggplot(df_ref, aes(x = DC1, y = DC2)) +
    geom_scattermore(aes(color = cell_type), pointsize = ptsize.ggscatermore) +
    # geom_point_rast(size=1) +
    theme_void() +
    coord_equal() +
    scale_color_manual(values = pal_hbca_siyuan)
# print(p0)

p10 <- p0 +
    new_scale_color() +
    geom_scattermore(
        data = df_pred,
        aes(x = DC1, y = DC2, color = archetype), pointsize = ptsize.ggscatermore
    ) +
    scale_color_manual(values = pal_ARC)
ggsave(
    filename = file.path(dir_res, "ref_diffusion_map.pdf"),
    plot = p0,
    width = 4, height = 4, useDingbats = FALSE
)
ggsave(
    filename = file.path(dir_res, "proj_diffusion_map.pdf"),
    plot = p10,
    width = 4, height = 4, useDingbats = FALSE
)

p1 <- ggplot(
    data = df_pred,
    aes(x = DC1, y = DC2, color = archetype)
) +
    geom_scattermore(pointsize = ptsize.ggscatermore) +
    coord_equal() +
    scale_color_manual(values = pal_ARC)
ggsave(
    filename = file.path(dir_res, "pred_diffusion_map.pdf"),
    plot = p1,
    width = 4, height = 4, useDingbats = FALSE
)
for (x in names(pal_ARC)) {
    p1x <- p0 +
        geom_scattermore(
            data = df_pred,
            aes(x = DC1, y = DC2), pointsize = ptsize.ggscatermore, color = "grey"
        ) +
        geom_scattermore(
            data = df_pred[df_pred$archetype %in% x, ],
            aes(x = DC1, y = DC2), pointsize = ptsize.ggscatermore, color = pal_ARC[x]
        ) +
        rremove("legend") + ggtitle(x)
    ggsave(
        filename = file.path(dir_res, paste0("pred_diffusion_map_", x, ".pdf")),
        plot = p1x,
        width = 4, height = 4, useDingbats = FALSE
    )
}

for (x in names(pal_pcr2)) {
    if (!x %in% df_pred$pCR_status) {
        next()
    }
    p1x <- p0 +
        geom_point_rast(
            data = df_pred,
            aes(x = DC1, y = DC2), size = ptsize.ggrastr, color = "grey"
        ) +
        geom_scattermore(
            data = df_pred[df_pred$pCR_status %in% x, ],
            aes(x = DC1, y = DC2), pointsize = ptsize.ggscatermore, color = pal_pcr2[x]
        ) +
        rremove("legend") + ggtitle(x)
    ggsave(
        filename = file.path(dir_res, paste0("pred_diffusion_map_", x, ".pdf")),
        plot = p1x,
        width = 4, height = 4, useDingbats = FALSE
    )
}

for (x in names(pal_ARC)) {
    for (y in names(pal_pcr2)) {
        idx <- df_pred$archetype %in% x & df_pred$pCR_status %in% y
        if (sum(idx) == 0) {
            next()
        }
        tmp <- df_pred[idx, ]
        tmp$pt_density <- get_density3(x = tmp[["DC1"]], y = tmp[["DC2"]])
        p2x <- p0 +
            geom_point_rast(
                data = df_pred[df_pred$archetype %in% x, ],
                aes(x = DC1, y = DC2), size = ptsize.ggrastr, color = "grey"
            ) +
            new_scale_color() +
            geom_scattermore(
                data = tmp,
                aes(x = DC1, y = DC2, color = pt_density), pointsize = ptsize.ggscatermore
            ) +
            scale_color_viridis_c(option = "E") +
            # scale_color_gradientn(colors = pal_viridis(10)) +
            rremove("legend") + ggtitle(paste0(x, " ", y))
        ggsave(
            filename = file.path(dir_res, paste0("pred_diffusion_map_", x, "_", y, ".pdf")),
            plot = p2x,
            width = 4, height = 4, useDingbats = FALSE
        )
    }
}


## is this patient specific?
pdf(file.path(dir_res, "pred_diffusion_map_patient.pdf"),
    onefile = TRUE,
    width = 4, height = 4, useDingbats = FALSE
)
tmp <- df_pred %>%
    dplyr::select(archetype, pCR_status, patient) %>%
    unique() %>%
    arrange(archetype, pCR_status) %>%
    dplyr::pull(patient)
dict_pat2arc <- df_meta %>%
    dplyr::select(patient, archetype) %>%
    unique() %>%
    deframe()
dict_pat2response <- df_meta %>%
    dplyr::select(patient, pCR_status) %>%
    unique() %>%
    deframe()
for (pat in tmp) {
    cat(match(pat, tmp), "/", length(tmp), "\n")
    ppat <- p0 +
        geom_point_rast(
            data = df_pred,
            aes(x = DC1, y = DC2), size = ptsize.ggrastr, color = "grey"
        ) + new_scale_color() +
        geom_scattermore(
            data = df_pred[df_pred$patient %in% pat, ],
            aes(x = DC1, y = DC2), pointsize = ptsize.ggscatermore * 2, color = "black"
        ) + rremove("legend") + ggtitle(sprintf(
            "%s %s %s (n=%d)", pat, dict_pat2arc[pat], dict_pat2response[pat],
            sum(df_pred$patient %in% pat)
        ))
    print(ppat)
}
dev.off()

#------------------ ~~~ Linear visualiation (pstime) ~~~ --------------------
cli_h1("Linear visualiation (pstime)")
library(patchwork)
library(rstatix)
library(ggpubr)
## use each DC component as a linear dimension to visualize the cells
dc <- "DC1"
for (dc in c(paste0("DC", 1:3))) {
    cli_alert_info("Visualizing {dc}")
    if (!dc %in% colnames(df_ref)) {
        cli_alert_warning("Skip {dc} as it is not in df_ref")
        next()
    }
    dc_vrange <- range(df_ref[[dc]])
    pl0 <- ggplot(df_ref, aes_string(x = dc, y = "cell_type")) +
        geom_violin(scale = "width", width = .9) +
        geom_boxplot(outlier.shape = NA, width = .2, aes(fill = cell_type)) +
        scale_fill_manual(values = pal_hbca_siyuan) +
        scale_x_continuous(limits = dc_vrange) +
        rremove("legend")
    pl1 <- ggplot(df_pred, aes_string(x = dc, y = "archetype")) +
        geom_violin(scale = "width", width = .9) +
        geom_boxplot(outlier.shape = NA, width = .2, aes(fill = archetype)) +
        scale_fill_manual(values = pal_ARC) +
        scale_x_continuous(limits = dc_vrange) +
        scale_y_discrete(limits = rev) +
        rremove("legend")
    pl <- pl0 + pl1 +
        plot_layout(ncol = 1, height = c(3, 4))
    # pl
    # pl <- ggarrange(pl0, pl1, ncol=1, align = 'v')
    ggsave(
        filename = file.path(dir_res, paste0("linear_visualization_", dc, ".pdf")),
        plot = pl,
        width = 3, height = 4, useDingbats = FALSE
    )
}

df_pred$archetype <- as.factor(df_pred$archetype)
df_pred$pCR_status <- fct_drop(df_pred$pCR_status)

for (dc in c(paste0("DC", 1:3))) {
    cli_alert_info("Visualizing {dc}")
    if (!dc %in% colnames(df_ref)) {
        cli_alert_warning("Skip {dc} as it is not in df_ref")
        next()
    }
    dc_vrange <- range(df_ref[[dc]])

    df_pred$value <- df_pred[[dc]]
    stat.test <- df_pred %>%
        group_by(archetype) %>%
        rstatix::wilcox_test(value ~ pCR_status, p.adjust.method = "fdr") %>% 
        adjust_pvalue(method = "BH") %>% 
        mutate(p.adj.signif = case_when(
            p.adj < 0.001 ~ "***",
            p.adj < 0.01 ~ "**",
            p.adj < 0.05 ~ "*",
            TRUE ~ "ns"
        ))
    stat.test$p.adj.signif[stat.test$p.adj.signif == "ns"] <- "" # remove ns
    stat.test$q_txt <- sprintf('q=%s', signif(stat.test$p.adj, digits = 3))
    stat.test$p_txt <- sprintf('p=%s', signif(stat.test$p, digits = 3))
    write.csv(
        stat.test,
        file = file.path(dir_res, paste0("linear_visualization_", dc, ".pCR_status.stat_test.csv")),
        row.names = FALSE
    )

    pl0 <- ggplot(df_ref, aes_string(x = dc, y = "cell_type")) +
        geom_violin(scale = "width", width = .9) +
        geom_boxplot(outlier.shape = NA, width = .2, aes(fill = cell_type)) +
        scale_fill_manual(values = pal_hbca_siyuan) +
        scale_x_continuous(limits = dc_vrange) +
        rremove("legend")
    pl1 <- ggplot(df_pred, aes_string(x = dc, y = "archetype")) +
        geom_violin(scale = "width", width = .8, aes(color = pCR_status), position = position_dodge(.8)) +
        geom_boxplot(outlier.shape = NA, width = .2, aes(fill = pCR_status), position = position_dodge(.8)) +
        geom_text(data=stat.test, aes(y=archetype, x=max(df_pred[[dc]]), label=q_txt), size=6/.pt) +
        scale_fill_manual(values = pal_pcr2) +
        scale_color_manual(values = pal_pcr2) +
        scale_y_discrete(limits = rev) +
        rremove("legend")
    ggsave(
        filename = file.path(dir_res, paste0("linear_visualization_", dc, ".pred_alone.pCR_status.pdf")),
        plot = pl1,
        width = 3, height = 3, useDingbats = FALSE
    )
    pl <- pl0 + (pl1 + scale_x_continuous(limits = dc_vrange)) +
        plot_layout(ncol = 1, height = c(3, 4 * 1.5))
    ggsave(
        filename = file.path(dir_res, paste0("linear_visualization_", dc, ".pCR_status.pdf")),
        plot = pl,
        width = 3, height = 5, useDingbats = FALSE
    )
    df_pred$value <- NULL
}








if (F) {
    # https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html
    # dm: DiffusionMap object
    ggplot(
        df_ref,
        aes(
            x = pseudotime_diffusionmap,
            y = cell_type2, colour = cell_type2
        )
    ) +
        geom_quasirandom(groupOnX = FALSE) +
        scale_color_manual(values = my_color) +
        theme_classic() +
        xlab("Diffusion map pseudotime (first diffusion map component)") +
        ylab("Timepoint") +
        ggtitle("Cells ordered by diffusion map pseudotime")
}
