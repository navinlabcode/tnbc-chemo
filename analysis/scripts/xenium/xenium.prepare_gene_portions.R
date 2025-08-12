suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Prepare 10 gene portions for xenium 5k or 5kPlus to facilitate
# evaluting label transfer (Seurat/Tangram)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: Oct 6, 2024
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
    library(ruok)
    library(Seurat)
    library(MASS)
    get_density2 <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95, fillna = 0) {
        # modified from http://slowkow.com/notes/ggplot2-color-by-density/
        # https://rdrr.io/github/GreenleafLab/ArchR/src/R/GgplotUtils.R
        x[is.na(x)] <- fillna
        y[is.na(y)] <- fillna
        df <- data.frame(x = x, y = y)
        dens <- MASS::kde2d(x = x, y = y, n = n)
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        df$density <- dens$z[ii]
        df$density[df$density > quantile(unique(df$density), densityMax)] <- quantile(unique(df$density), densityMax) # make sure the higher end doesnt bias colors
        if (!is.null(sample)) {
            df <- df[sample(nrow(df), min(sample, nrow(df))), ]
        }
        return(df)
    }
})
cmdargs <- commandArgs(trailingOnly = TRUE)

if (length(cmdargs) > 0) {
    f_ref <- cmdargs[1]
    f_sp <- cmdargs[2]
    xenium_type <- cmdargs[3] # xenium5k xenium5kPlus
    dir_res <- file.path(dirname(f_ref), paste0("gene_portions_", xenium_type))
} else {
	## Input to the Seurat object of scRNA
    f_ref <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/downto_1000_by_cell_state_paper/ready.sr3.rds"
    f_ref <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/T/downto_1000_by_cell_state_paper/ready.sr3.rds"

    ## xenium5k gene libraray
    xenium_type <- "xenium5k"
    f_sp <- "/volumes/spatial/xenium/20240809__195707__5kpanel_BCMDCIS41T2-1_ECIS06T2_ART10pretx_ART312pretx_ART23pretx_ART304preTX_ART305pretx_ART311pretx_BCMDCIS05T1_ECIS23T2_BCMDCIS42T_R2_ECIS05T-S1_5K_08092024/output-XETG00093__0037629__ART10PreTx_5K__20240809__200032/cell_feature_matrix/features.tsv.gz"
    dir_res <- file.path(dirname(f_ref), "gene_portions_xenium5k")

    ## xenium5kPlus gene library
    # f_sp <- '/volumes/spatial/xenium/20240911__191532__5k_ARTADD_ON_panel_ART18pretx_ART31pretx_ART43pretx_ART65pretx_ART40pretx_ESOp4pre_mid_ESOp36pre_mid_ESOp24pre_mid__BCMDCIS102T3_BCMDCIS31T_ESOp37_pre_mid_ESOp8_pre_mid_ESOp37pre_mid_ESOp44_pre_mid_ESOp23_pre_mid_OCT_09112024/output-XETG00093__0041100__ART18pretx__20240911__192016/cell_feature_matrix/features.tsv.gz'
    # dir_res <- file.path(dirname(f_ref), 'gene_portions_xenium5kPlus')
}

#------------------ ~~~ Rock ~~~ --------------------
dir_create(dir_res)

obj <- read_rds(f_ref)
print(obj)

DefaultAssay(obj) <- "RNA"
ref_genes_symbol <- rownames(obj)
str(ref_genes_symbol)

sp_genes <- read_tsv(f_sp, col_names = F)
head(sp_genes)
dim(sp_genes)

## xenium genes (Use 'Gene Expression'. Others include the positive/negative control probes)
table(sp_genes$X3)
sp_genes <- sp_genes[sp_genes$X3 == "Gene Expression", ]
sp_genes_symbol <- as.character(sp_genes$X2)
str(sp_genes_symbol)

if (!all(sp_genes_symbol %in% ref_genes_symbol)) {
    cli_alert_warning("not all spatial genes are matched with the single-cell reference; probably due to version difference.")
    unmatched_sp_genes <- sp_genes_symbol[!sp_genes_symbol %in% ref_genes_symbol]
    str(unmatched_sp_genes)
    write_lines(
        unmatched_sp_genes,
        file.path(dir_res, "spatial_gene_symbols_unmatching_sc.txt")
    )
}

genes_shared <- intersect(sp_genes_symbol, ref_genes_symbol)
str(genes_shared)
write_lines(
    genes_shared,
    file.path(dir_res, "spatial_gene_symbols_matching_sc.txt")
)

#------------------ ~~~ Check HVFinfo ~~~ --------------------
##
## Check if the spatial gene panels have the similar stats pattern as the
## original entire gene sets.
##
DefaultAssay(obj) <- "RNA"
df_ref_hv <- HVFInfo(obj)
dim(df_ref_hv)
all(ref_genes_symbol %in% rownames(df_ref_hv))
# print(VariableFeaturePlot(obj))
df_ref_hv$is_in_spatial <- as.factor(ifelse(rownames(df_ref_hv) %in% genes_shared,
    "spatial_yes", "spatial_no"
))
table(df_ref_hv$is_in_spatial)

# library(scattermore)
library(ggrastr)
head(df_ref_hv)
df_ref_hv <- rownames_to_column(df_ref_hv, "gene")
## add scaterplot density
df_ref_hv <- lapply(
    split(df_ref_hv, df_ref_hv$is_in_spatial),
    function(df) {
        tmp <- get_density2(x = log10(df$mean + 1), y = df$variance.standardized, n = 200)
        df$scater_density <- tmp$density
        df$scater_density <- df$scater_density / max(df$scater_density)
        return(df)
    }
) %>% do.call("rbind", .)
rownames(df_ref_hv) <- df_ref_hv$gene

p <- df_ref_hv %>%
    dplyr::arrange(scater_density) %>%
    ggplot(., aes_string(
        x = "mean", y = "variance.standardized",
        fill = "scater_density"
    )) +
    geom_point_rast(pch = 21, color = "lightgrey", stroke = .1) +
    # facet_grid(cols = vars(is_in_spatial)) +
    facet_wrap(~is_in_spatial, nrow = 1) +
    scale_x_log10() +
    labs(x = "log10(mean expression)", y = "standardized variance") +
    scale_fill_viridis_c(option = "E")
# print(p)
ggsave(file.path(dir_res, "scatter.HVGplot.by_in_spatial_panel.pdf"), p, width = 6, height = 3, useDingbats = F)

## Label some known genes to have some understanding of the scatter plot.
marker_manual <- c(
    "EPCAM", "TACSTD2",
    "PTPRC", "CD68", "SPI1",
    "CD3G", "CD79A",
    "LUM", "DCN", "PECAM1", "VWF", "RGS5", "STEAP4"
)
marker_data <- df_ref_hv %>%
    dplyr::group_by(is_in_spatial) %>%
    dplyr::top_n(wt = variance.standardized, n = 3) %>%
    dplyr::ungroup() %>%
    dplyr::pull(gene)

marker_tolabel <- unique(c(marker_manual, marker_data))

df_marker_tolabel <- df_ref_hv[marker_tolabel, ]
df_marker_tolabel <- df_marker_tolabel[df_marker_tolabel$mean != 0, ]
marker_tolabel <- df_marker_tolabel$gene

head(df_marker_tolabel)
library(ggrepel)

p <- df_ref_hv %>%
    dplyr::arrange(scater_density) %>%
    ggplot(., aes_string(
        x = "mean", y = "variance.standardized",
        fill = "scater_density"
    )) +
    geom_point_rast(pch = 21, color = "lightgrey", stroke = .1) +
    geom_point(
        data = df_marker_tolabel, aes_string(
            x = "mean", y = "variance.standardized"
        ),
        pch = 21, fill = NA, color = "red", stroke = .5
    ) +
    ggrepel::geom_text_repel(
        data = df_marker_tolabel, aes_string(
            x = "mean", y = "variance.standardized", label = "gene"
        ),
        min.segment.length = unit(0, "lines"), color = "red",
        size = 4 / .pt, max.overlaps = 100
    ) +
    facet_grid(cols = vars(is_in_spatial)) +
    scale_x_log10() +
    labs(x = "log10(mean expression)", y = "standardized variance") +
    scale_fill_viridis_c(option = "E")
# print(p)
ggsave(file.path(dir_res, "scatter.HVGplot.by_in_spatial_panel_label.pdf"),
    p,
    width = 6, height = 3, useDingbats = F
)

## in case of future needs to follow up with some new analysis, load this object
write_rds(df_ref_hv, file.path(dir_res, "rdata.df_ref_hv.df.rds"))


gene_nondetected <- df_ref_hv$gene[df_ref_hv$mean == 0]
gene_nonvariable <- df_ref_hv$gene[df_ref_hv$variance.standardized == 0]
str(genes_shared)
str(gene_nonvariable)

write_lines(
    intersect(gene_nondetected, genes_shared),
    file.path(dir_res, "spatial_gene_symbols_matching_sc.nondetected_in_sc.txt")
)
write_lines(
    intersect(gene_nonvariable, genes_shared),
    file.path(dir_res, "spatial_gene_symbols_matching_sc.nonvariable_in_sc.txt")
)

genes_use <- setdiff(genes_shared, gene_nondetected)
genes_use <- setdiff(genes_use, gene_nonvariable)
str(genes_use)

library(UpSetR)
upset_list <- list(
    SC = ref_genes_symbol,
    SP = sp_genes_symbol,
    Not_detected_SC = gene_nondetected,
    Not_variable_SC = gene_nonvariable,
    InUse = genes_use
)

pdf(file.path(dir_res, "upsetplot.genes_filtering_progress.pdf"), width = 4, height = 4, useDingbats = F)
upset(fromList(upset_list),
    sets = rev(names(upset_list)), mb.ratio = c(0.5, 0.5),
    keep.order = TRUE
)
dev.off()


df_ref_hv_sp <- df_ref_hv[genes_use, ]
dim(df_ref_hv_sp)

#------------------ ~~~ Split genes ~~~ --------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Ensure each gene portion has the similar stats property of the
## entire spatial gene panel.

## All genes are split into 10 bins based on log10(mean).
## Genes each bin are split into 10 portions.
## Combine genes of each portion.

## random sampling may not guarantee this requirement.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(df_ref_hv_sp)
range(log10(df_ref_hv_sp$mean))
hist(log10(df_ref_hv_sp$mean))
df_ref_hv_sp$bins <- cut(log10(df_ref_hv_sp$mean),
    breaks = 10, right = FALSE
)
df_ref_hv_sp$bins <- as.factor(as.numeric(df_ref_hv_sp$bins))

N_PORTIONS <- 10
set.seed(1026)
df_ref_hv_sp <- lapply(
    split(df_ref_hv_sp, df_ref_hv_sp$bins),
    function(df) {
        print(nrow(df))
        i_shuffled <- sample(seq_len(nrow(df)), size = nrow(df))
        i_portion <- i_shuffled %% N_PORTIONS ## 1,2,..,8,9,0
        i_portion[i_portion == 0] <- N_PORTIONS
        df$portion <- i_portion
        return(df)
    }
) %>% do.call("rbind", .)
rownames(df_ref_hv_sp) <- df_ref_hv_sp$gene
stopifnot(identical(sort(genes_use), sort(rownames(df_ref_hv_sp))))
df_ref_hv_sp <- df_ref_hv_sp[genes_use, ]


tmp <- get_density2(
    x = log10(df_ref_hv_sp$mean),
    y = df_ref_hv_sp$variance.standardized, n = 100
)

df_ref_hv_sp$scater_density <- tmp$density
df_ref_hv_sp$scater_density <- df_ref_hv_sp$scater_density / max(df_ref_hv_sp$scater_density, na.rm = T)


p <- df_ref_hv_sp %>%
    dplyr::arrange(scater_density) %>%
    ggplot(., aes_string(
        x = "mean", y = "variance.standardized",
        fill = "scater_density"
    )) +
    geom_point_rast(pch = 21, color = "lightgrey", stroke = .1) +
    # facet_wrap(~is_in_spatial, nrow = 1) +
    scale_x_log10() +
    labs(x = "log10(mean expression)", y = "standardized variance") +
    scale_fill_viridis_c(option = "E") +
    theme(aspect.ratio = 1)
# print(p)
ggsave(file.path(dir_res, "scatter.HVGplot.genes_in_use.pdf"), p, width = 4, height = 3, useDingbats = F)


table(df_ref_hv_sp$portion)

df_ref_hv_sp <- lapply(
    split(df_ref_hv_sp, df_ref_hv_sp$portion),
    function(df) {
        tmp <- get_density2(
            x = log10(df$mean),
            y = df$variance.standardized, n = 100
        )
        df$scater_density_portion <- tmp$density
        df$scater_density_portion <- df$scater_density_portion / max(df$scater_density_portion)
        return(df)
    }
) %>% do.call("rbind", .)
rownames(df_ref_hv_sp) <- df_ref_hv_sp$gene

p <- df_ref_hv_sp %>%
    dplyr::group_by(portion) %>%
    dplyr::arrange(scater_density, .by_group = TRUE) %>%
    ggplot(., aes_string(
        x = "mean", y = "variance.standardized",
        fill = "scater_density"
    )) +
    geom_point_rast(pch = 21, color = "lightgrey", stroke = .1) +
    facet_wrap(~portion, nrow = 2) +
    scale_x_log10() +
    labs(x = "log10(mean expression)", y = "standardized variance") +
    scale_fill_viridis_c(option = "E") +
    theme(aspect.ratio = 1)
# print(p)
ggsave(file.path(dir_res, "scatter.HVGplot.genes_in_use.by_portions.pdf"),
    p + rremove("legend"),
    width = 5 * 2, height = 2 * 2, useDingbats = F
)

## in case of future needs to follow up with some new analysis, load this object
write_rds(df_ref_hv_sp, file.path(dir_res, "rdata.df_ref_hv_sp.df.rds"))
write_csv(df_ref_hv_sp, file.path(dir_res, "rdata.df_ref_hv_sp.df.csv"))

head(df_ref_hv_sp)

for (df in split(df_ref_hv_sp, df_ref_hv_sp$portion)) {
    portion_name <- unique(df$portion)

    write_rds(
        as.character(df$gene),
        file.path(dir_res, sprintf("genes_of_portion_%s.rds", portion_name))
    )
    write_lines(
        as.character(df$gene),
        file.path(dir_res, sprintf("genes_of_portion_%s.txt", portion_name))
    )
    cli_alert_success(c("Portion-", portion_name))
    str(as.character(df$gene))
}
rm(df)

cat("[done R xenium.prepare_gene_portions.R]")
