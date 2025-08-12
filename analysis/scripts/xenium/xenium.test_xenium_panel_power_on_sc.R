suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Test whether the xenium panel is powerful enough to
# recapitulate the cell states (`cat_transfer`) in each cell compartment.
#
# Create a new assay XNA in the object by using the shared
# genes between the original assay (`reference.assay`) and the
# xenium panel (`xenium_panel_type`, either the 5k or the 5k+addons)
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: text
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
    my_scatter_themevoid <- theme_pubr(base_size = 6, legend = "right") %+replace% theme(
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = rel(1)),
        axis.line = element_blank()
    )
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
    # library(RhpcBLASctl)
    # RhpcBLASctl::blas_set_num_threads(8)
    # library(future)
    # plan("multisession", workers = 8)
    # plan()
    # options(future.globals.maxSize = 20 * 1024^3) # for 20 Gb RAM
})
cmdargs <- commandArgs(trailingOnly = TRUE)

if (length(cmdargs) > 0) {
    f_ref <- cmdargs[1]
    reference.assay <- cmdargs[2]
    cat_transfer <- cmdargs[3]
    xenium_panel_type <- cmdargs[4]
} else {
    ## ref
    # f_ref           <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/T/downto_1000_by_cell_state_paper/ready.sr3.rds'
    # cat_transfer    <- 'cell_state_paper'
    # reference.assay <-  'RNA'
    # xenium_panel_type <- 'xenium5k'

    f_ref <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Fibro/downto_1000_by_cell_state_paper/ready.sr3.rds"
    cat_transfer <- "cell_state_paper"
    reference.assay <- "RNA"
    xenium_panel_type <- "xenium5k"

    # obj <- read_rds("/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/objects_split_into_celltype/Fibro/downto_1000_by_cell_state_paper/test_power_xenium5k/ready.seurat.rds")
}


message(f_ref)
message(reference.assay)
message(cat_transfer)
message(xenium_panel_type)

cell_type_compartment <- basename(dirname(dirname(f_ref)))
message(c(cell_type_compartment, " cells"))


f_sp <- switch(xenium_panel_type,
    "xenium5k" = "/volumes/spatial/xenium/20240809__195707__5kpanel_BCMDCIS41T2-1_ECIS06T2_ART10pretx_ART312pretx_ART23pretx_ART304preTX_ART305pretx_ART311pretx_BCMDCIS05T1_ECIS23T2_BCMDCIS42T_R2_ECIS05T-S1_5K_08092024/output-XETG00093__0037629__ART10PreTx_5K__20240809__200032/cell_feature_matrix/features.tsv.gz",
    "xenium5kPlus" = "/volumes/spatial/xenium/20240911__191532__5k_ARTADD_ON_panel_ART18pretx_ART31pretx_ART43pretx_ART65pretx_ART40pretx_ESOp4pre_mid_ESOp36pre_mid_ESOp24pre_mid__BCMDCIS102T3_BCMDCIS31T_ESOp37_pre_mid_ESOp8_pre_mid_ESOp37pre_mid_ESOp44_pre_mid_ESOp23_pre_mid_OCT_09112024/output-XETG00093__0041100__ART18pretx__20240911__192016/cell_feature_matrix/features.tsv.gz"
)

dir_res <- file.path(dirname(f_ref), sprintf("test_power_%s", xenium_panel_type))
fs::dir_create(dir_res)

if (TRUE) {
    ## read xenium genes
    sp_genes <- read_tsv(f_sp, col_names = F)
    head(sp_genes)
    dim(sp_genes)
    ## xenium genes (Use 'Gene Expression'. Others include the positive/negative control probes)
    table(sp_genes$X3)
    sp_genes <- sp_genes[sp_genes$X3 == "Gene Expression", ]
    genes_xenium <- as.character(sp_genes$X2)
    str(genes_xenium)
}

std_ident_levels <- read_rds(
    file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate",
        "pat102/atlas",
        "idents_levels_cellstates.rds"
    )
)

obj <- read_rds(f_ref)
print(obj)

obj@meta.data[, cat_transfer] <- standardize_factor(obj@meta.data[, cat_transfer], std_ident_levels)
Idents(obj) <- cat_transfer
table(Idents(obj), useNA = "ifany")
pal_ident <- init_pal_d(levels(Idents(obj)))

ggsave(
    file.path(dir_res, sprintf("legend.%s.pdf", cat_transfer)),
    pal_to_ggplot(pal_ident, pal_name = cat_transfer),
    width = 3, height = 3, useDingbats = F
)


p0 <- DimPlot(obj,
    reduction = "umap", label = T, shuffle = T,
    raster = T, raster.dpi = c(1024, 1024), pt.size = 3
) +
    my_scatter_themevoid +
    scale_color_manual(values = pal_ident)
print(p0)

pdf(file.path(dir_res, sprintf("dimplot.%s.%s.pdf", "umap_atlas", cat_transfer)),
    width = 6, height = 6, useDingbats = F
)
print(p0 + rremove("legend"))
dev.off()
#------------------ ~~~ Shard genes ~~~ --------------------

genes_sc <- rownames(obj)
str(genes_sc)
genes_shared <- intersect(genes_xenium, genes_sc)
setdiff(genes_xenium, genes_sc)

str(genes_shared)
n_genes_shared <- length(genes_shared)
#------------------ ~~~ XNA: xenium&RNA shared ~~~ --------------------

assay_xna <- obj[["RNA"]]
assay_xna <- subset(assay_xna, features = genes_shared)
assay_xna <- CreateSeuratObject(counts = assay_xna@counts, assay = "XNA")
obj[["XNA"]] <- assay_xna[["XNA"]]
obj <- AddMetaData(obj, metadata = assay_xna@meta.data[, c("nCount_XNA", "nFeature_XNA")])

DefaultAssay(obj) <- "XNA"
obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3000)

if (FALSE) {
    # probabaly not recommend because using 3000 hvg is similar to the original single-cell.
    is_non_empty_features <- rowSums(GetAssayData(obj, layer = "counts", assay = "XNA")) != 0
    print(table(is_non_empty_features))
    non_empty_features <- rownames(GetAssayData(obj, layer = "counts", assay = "XNA"))[is_non_empty_features]
    VariableFeatures(obj) <- non_empty_features
}

obj <- ScaleData(obj, vars.to.regress = "nCount_XNA")
obj <- RunPCA(obj,
    npcs = 100,
    reduction.name = "pca.xna", reduction.key = "PCXNA"
)
obj <- RunUMAP(obj,
    dims = 1:50,
    reduction.name = "umap.xna", reduction.key = "UMAPXNA",
    reduction = "pca.xna", return.model = TRUE
)

pp <- DimPlot(obj,
    reduction = "umap.xna", label = T, shuffle = T,
    raster = T, raster.dpi = c(1024, 1024), pt.size = 3
) +
    my_scatter_themevoid +
    scale_color_manual(values = pal_ident)
print(pp)

pdf(file.path(dir_res, sprintf("dimplot.%s.%s.pdf", "umap_xna", cat_transfer)),
    width = 6, height = 6, useDingbats = F
)
print(pp + rremove("legend"))
dev.off()

xna_hvg <- VariableFeatures(obj)
write_lines(xna_hvg, file.path(dir_res, "xna.hvg_names.txt"))
write_rds(xna_hvg, file.path(dir_res, "xna.hvg_names.rds"))

#------------------ ~~~ Original RNA ~~~ --------------------
DefaultAssay(obj) <- "RNA"
obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData(vars.to.regress = "nCount_RNA")
obj <- RunPCA(obj, npcs = 100)
obj <- RunUMAP(obj, dims = 1:50, return.model = TRUE)
pp0 <- DimPlot(obj,
    reduction = "umap", label = T, shuffle = T,
    raster = T, raster.dpi = c(1024, 1024), pt.size = 3
) +
    my_scatter_themevoid +
    scale_color_manual(values = pal_ident)
print(pp0)

pdf(file.path(dir_res, sprintf("dimplot.%s.%s.pdf", "umap_perse", cat_transfer)),
    width = 6, height = 6, useDingbats = F
)
print(pp0 + rremove("legend"))
dev.off()

#------------------ ~~~ Export ~~~ --------------------
DefaultAssay(obj) <- "XNA"
stopifnot(identical(VariableFeatures(obj), xna_hvg)) # seurat v5 stores hvg for different assays
print(obj)
# An object of class Seurat
# 38502 features across 14000 samples within 2 assays
# Active assay: XNA (4964 features, 3000 variable features)
# 1 other assay present: RNA
# 4 dimensional reductions calculated: pca, umap, pca.xna, umap.xna
print(table(Idents(obj)))
write_seurat(obj, dir_res, "ready")

cat("[done]")
timestamp()
