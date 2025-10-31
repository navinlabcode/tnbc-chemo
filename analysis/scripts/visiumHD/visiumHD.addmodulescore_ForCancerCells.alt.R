suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Add module scores to cancer cells only.
#
# cell identity is stored in a separate file.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-08-17 / 2024-12-09 / 2025-02-26 / 2025-05-14
# Based on 'visiumHD.addmodulescore_ForCancerSubset.R'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
" -> doc_help
timestamp()
suppressPackageStartupMessages({
    library(readr)
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    theme_set(theme_pubr(base_size = 8, legend = "right") %+replace% theme(axis.ticks.length = unit(0.1, "inch")))
    library(cli)
    library(tictoc)
    library(glue)
    library(scales)
    library(tools)
    library(ruok)
    library(Seurat)
    # library(arrow)
    library(fs)
    library(Signac)
    library(clustree)
    library(UCell)
    library(AUCell)
    library(BiocParallel)
    my_scatter_themevoid <- theme_pubr(base_size = 8, legend = "right") %+replace% theme(
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = rel(1)),
        axis.line = element_blank()
    )
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("/volumes/USR1/yyan/apps/bonnie/R/io_gsea_formats.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/uti.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
})
options <- commandArgs(trailingOnly = TRUE)

# Example
# f_in <- "/volumes/USR1/yyan/project/tnbc_visium_hd/data/ART258/binsize_32/ready.seurat.rds"
# f_df <- "/volumes/USR1/yyan/project/tnbc_visium_hd/data/ART258/binsize_32/finalize_cancer_cells/celltype_final.rds"
# cellmeta_suit_name <- "css_celltypes"
# ident_use <- "celltypes"
# cancer_ident <- "Tumor"
# gene_suit_name <- c(
#     "ARC", "ARC_updated",
#     "ARC_updated_genes_xenium5k", "ARC_updated_genes_xenium5kPlus",
#     "MP"
# )
# gene_suit_name <- "ARC_updated_genes_xenium5kPlus"

if (length(options) > 0) {
    f_in <- options[1]
    f_df <- options[2]
    cellmeta_suit_name <- options[3]
    ident_use <- options[4]
    cancer_ident <- options[5]
    gene_suit_name <- options[6]
}

cli_alert_info(gene_suit_name)

message(f_in)
xmo <- read_rds(f_in)
print(xmo)

df_in <- read_rds(f_df)
print(head(df_in))

if ("Barcode" %in% colnames(df_in)) {
    df_in <- as.data.frame(df_in)
    rownames(df_in) <- df_in$Barcode
}

xmo <- AddMetaData(xmo, df_in)
stopifnot(ident_use %in% colnames(df_in))
stopifnot(ident_use %in% colnames(xmo@meta.data))
stopifnot(cancer_ident %in% xmo@meta.data[, ident_use])

xmo$celltype <- xmo[[]][, ident_use]
ident_use <- "celltype"

try(FetchData(xmo, vars = "EPCAM") %>% deframe() %>% range() %>% print())

# method_opts <- c("module", "AMSREV", "ucell", "aucell")
method_opts <- c("module", "aucell")
cli_li(method_opts)
force_do_module <- TRUE
force_do_amsrev <- TRUE
force_do_ucell <- TRUE
force_do_aucell <- TRUE

#------------------ ~~~ lib: gene lists ~~~ --------------------
if (T) {
    #------ Archetype old based on NMF only ------
    arc_genes_old <- read_rds("~/project/tnbc_pre_atlas//rds_rna-integrate/pat102/lv01.aneuploidy_tri_type.aneuploid.pure5/psbulk/fastnmf/rank4/deliver.nmf_markers.top.rds")
    names(arc_genes_old) <- gsub(pattern = "fNMF", replacement = "ARC", names(arc_genes_old))
    str(arc_genes_old)

    #------ Archetype updated: NMF + DEG ------
    arc_genes <- read_rds("/volumes/USR1/yyan/project/tnbc_pre_atlas/summary_fig/archetype_psbulkDEG_MAplot/fnmf_marker_sig.list.rds")
    names(arc_genes) <- gsub(pattern = "fNMF", replacement = "ARC", names(arc_genes))
    str(arc_genes)

    #------ Archetype update x xenium5k ------
    arc_genes_x5k <- read_rds("/volumes/USR1/yyan/project/tnbc_xenium/lib/xenium_probes/xenium5k_genes.rds")
    arc_genes_x5k <- lapply(arc_genes, function(v) {
        intersect(v, arc_genes_x5k)
    })
    str(arc_genes_x5k)
    #------ Archetype update x xenium5kPlus ------
    arc_genes_x5kPlus <- read_rds("/volumes/USR1/yyan/project/tnbc_xenium/lib/xenium_probes/xenium5kPlus_genes.rds")
    arc_genes_x5kPlus <- lapply(arc_genes, function(v) {
        intersect(v, arc_genes_x5kPlus)
    })
    str(arc_genes_x5kPlus)

    #------ Archetype DEGsLv30 x xenium5k ------
    arc_DEGsLv30_x5k <- read_rds("/volumes/USR1/yyan/project/tnbc_pre_atlas/summary_fig/archetype_psbulkDEG_MAplot/archetype_DEGsLv30_overlaping_xenium_probes/arc_genes_xenium5k.rds")
    str(arc_DEGsLv30_x5k)
    #------ Archetype DEGsLv30 x xenium5kPlus ------
    arc_DEGsLv30_x5kPlus <- read_rds("/volumes/USR1/yyan/project/tnbc_pre_atlas/summary_fig/archetype_psbulkDEG_MAplot/archetype_DEGsLv30_overlaping_xenium_probes/arc_genes_xenium5kPlus.rds")
    str(arc_DEGsLv30_x5kPlus)
    #------ cancer metaprogram ------
    metaprogram_genes <- read_rds(
        file.path(
            "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate",
            "pat102/lv01.aneuploidy_tri_type.aneuploid.pure5",
            "metamodule_fnmf",
            "MM_alt_clean_byscore_wardD2",
            "deliver.mm_markers.rds"
        )
    )
    metaprogram_genes <- metaprogram_genes[c(1, 4:13)]
    str(metaprogram_genes)
}


#------------------ ~~~ setup output dir ~~~ --------------------
dir_res <- file.path(
    dirname(f_in),
    sprintf("%s.Tumor_cells_only_modulescore_%s", cellmeta_suit_name, gene_suit_name)
)
# dir_res <- file.path(
#     dirname(f_in),
#     sprintf("Tumor_cells_only_modulescore_%s", gene_suit_name)
# )

message(dir_res)
fs::dir_create(dir_res)


gene_set_use <- switch(gene_suit_name,
    ARC = arc_genes_old,
    ARC_updated = arc_genes,
    ARC_updated_genes_xenium5k = arc_genes_x5k,
    ARC_updated_genes_xenium5kPlus = arc_genes_x5kPlus,
    ARC_DEGs_xenium5k = arc_DEGs_x5k,
    ARC_DEGs_xenium5kPlus = arc_DEGs_x5kPlus,
    ARC_DEGsLv30_xenium5k = arc_DEGsLv30_x5k,
    ARC_DEGsLv30_xenium5kPlus = arc_DEGsLv30_x5kPlus,
    MP = metaprogram_genes
)


gene_set_use <- lapply(gene_set_use, function(v) {
    intersect(v, rownames(xmo))
})
str(gene_set_use)

gene_set_use_size <- sapply(gene_set_use, length)

if (sum(gene_set_use_size == 0) > 0) {
    cli_alert_danger("These genesets have no genes found in the object: ")
    cli_ol(names(gene_set_use_size)[gene_set_use_size == 0])
    gene_set_use <- gene_set_use[gene_set_use_size > 0]
}


write_rds(gene_set_use, file.path(dir_res, "gene_set.list.rds"))
write_lines(unlist(gene_set_use), file.path(dir_res, "gene_set.txt"))
write_gmx2(gene_set_use, file.path(dir_res, "gene_set.gmx.csv"))
capture.output(str(gene_set_use), file = file.path(dir_res, "gene_set.glimpse.txt"))

table(xmo$celltype, useNA = "ifany") %>%
    c() %>%
    enframe() %>%
    write_csv(file.path(dir_res, "celltype_counts.csv"))

#------ only calc on cancer cells ------

cli_alert_info(c("cancer ident is: ", ident_use, "==", cancer_ident))

Idents(xmo) <- ident_use
idx <- c(Idents(xmo) == cancer_ident)
table(idx)
if (!all(idx)) {
    message("subsetting cell...")
    cancer_cnames <- Cells(xmo)[idx]
    str(cancer_cnames)
    xmo_sub <- subset(xmo, cells = cancer_cnames)
} else {
    xmo_sub <- xmo
}

cli_alert_info(c("Computing ", comma(ncol(xmo_sub)), " cells in total. "))
#------------------ ~~~ AddModuleScore ~~~ --------------------
cli_h1("AddModuleScore")
#------ helper functions ------
if (T) {
    AMS_rev <- function(object, features, pool = NULL, nbin = 24, ctrl = 100,
                        k = FALSE, assay = NULL, name = "Cluster", seed = 1, search = FALSE,
                        ...) {
        if (!is.null(x = seed)) {
            set.seed(seed = seed)
        }
        assay.old <- DefaultAssay(object = object)
        assay <- assay %||% assay.old
        DefaultAssay(object = object) <- assay
        assay.data <- GetAssayData(object = object)
        features.old <- features
        if (k) {
            .NotYetUsed(arg = "k")
            features <- list()
            for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
                features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster ==
                    i))
            }
            cluster.length <- length(x = features)
        } else {
            if (is.null(x = features)) {
                stop("Missing input feature list")
            }
            features <- lapply(X = features, FUN = function(x) {
                missing.features <- setdiff(x = x, y = rownames(x = object))
                if (length(x = missing.features) > 0) {
                    warning("The following features are not present in the object: ",
                        paste(missing.features, collapse = ", "),
                        ifelse(test = search, yes = ", attempting to find updated synonyms",
                            no = ", not searching for symbol synonyms"
                        ),
                        call. = FALSE, immediate. = TRUE
                    )
                    if (search) {
                        tryCatch(expr = {
                            updated.features <- UpdateSymbolList(
                                symbols = missing.features,
                                ...
                            )
                            names(x = updated.features) <- missing.features
                            for (miss in names(x = updated.features)) {
                                index <- which(x == miss)
                                x[index] <- updated.features[miss]
                            }
                        }, error = function(...) {
                            warning("Could not reach HGNC's gene names database",
                                call. = FALSE, immediate. = TRUE
                            )
                        })
                        missing.features <- setdiff(x = x, y = rownames(x = object))
                        if (length(x = missing.features) > 0) {
                            warning("The following features are still not present in the object: ",
                                paste(missing.features, collapse = ", "),
                                call. = FALSE, immediate. = TRUE
                            )
                        }
                    }
                }
                return(intersect(x = x, y = rownames(x = object)))
            })
            cluster.length <- length(x = features)
        }
        if (!all(Seurat:::LengthCheck(values = features))) {
            warning(paste(
                "Could not find enough features in the object from the following feature lists:",
                paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
                "Attempting to match case..."
            ))
            features <- lapply(
                X = features.old, FUN = CaseMatch,
                match = rownames(x = object)
            )
        }
        if (!all(Seurat:::LengthCheck(values = features))) {
            stop(paste(
                "The following feature lists do not have enough features present in the object:",
                paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
                "exiting..."
            ))
        }
        pool <- pool %||% rownames(x = object)
        data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
        data.avg <- data.avg[order(data.avg)]
        data.cut <- cut_number(
            x = data.avg + rnorm(n = length(data.avg)) / 1e+30,
            n = nbin, labels = FALSE, right = FALSE
        )
        names(x = data.cut) <- names(x = data.avg)
        ctrl.use <- vector(mode = "list", length = cluster.length)
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            for (j in 1:length(x = features.use)) {
                ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut ==
                    data.cut[features.use[j]])], size = ctrl, replace = TRUE)))
            }
        }
        ctrl.use <- lapply(X = ctrl.use, FUN = unique)
        ctrl.scores <- matrix(
            data = numeric(length = 1L), nrow = length(x = ctrl.use),
            ncol = ncol(x = object)
        )
        for (i in 1:length(ctrl.use)) {
            features.use <- ctrl.use[[i]]
            ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
        }
        features.scores <- matrix(
            data = numeric(length = 1L), nrow = cluster.length,
            ncol = ncol(x = object)
        )
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            data.use <- assay.data[features.use, , drop = FALSE]
            features.scores[i, ] <- Matrix::colMeans(x = data.use)
        }
        features.scores.use <- features.scores - ctrl.scores
        rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
        features.scores.use <- as.data.frame(x = t(x = features.scores.use))
        rownames(x = features.scores.use) <- colnames(x = object)
        object[[colnames(x = features.scores.use)]] <- features.scores.use
        Seurat:::CheckGC()
        DefaultAssay(object = object) <- assay.old
        return(object)
    }
}

for (module_method in method_opts) {
    cli_h2(module_method)
    message(c("module method = ", module_method))
    f_module <- file.path(dir_res, sprintf("%sscore.df.rds", module_method))
    module_score_prefix <- paste0(module_method, "_score")

    if (file.exists(f_module)) {
        # if (F) {
        cat("Load the existing results... ")

        if (any(str_detect(colnames(xmo@meta.data), module_score_prefix))) {
            cat("the object had results so overwrite...")
        }

        df_module <- read_rds(f_module)
        stopifnot(identical(rownames(df_module), Cells(xmo)))
        xmo <- AddMetaData(xmo, df_module)
    } else {
        cat("Computing on ", cancer_ident, " cells only ...")
        run_module_func <- switch(module_method,
            module = AddModuleScore,
            AMSREV = AMS_rev,
            ucell = AddModuleScore_UCell
        )
        module_score_prefix <- paste0(module_method, "_score")
        name_module_set <- paste0(module_score_prefix, "_", names(gene_set_use))

        if (module_method %in% c("module", "AMSREV")) {
            #------ Tirosh or AMSREV ------
            xmo_sub <- run_module_func(xmo_sub, gene_set_use,
                name = module_score_prefix
            )
            for (i in seq_along(names(gene_set_use))) {
                i_old <- paste0(module_score_prefix, i)
                i_new <- name_module_set[i]
                xmo_sub@meta.data[[i_new]] <- xmo_sub@meta.data[[i_old]]
                xmo_sub@meta.data[[i_old]] <- NULL
            }
        } else if (module_method == "ucell") {
            #------ UCELL ------
            xmo_sub <- run_module_func(xmo_sub, gene_set_use, name = "_UCell")
            for (i in seq_along(names(gene_set_use))) {
                i_old <- paste0(names(gene_set_use)[i], "_UCell")
                i_new <- name_module_set[i]
                xmo_sub@meta.data[[i_new]] <- xmo_sub@meta.data[[i_old]]
                xmo_sub@meta.data[[i_old]] <- NULL
            }
        } else if (module_method == "aucell") {
            #------ AUCELL ------
            cat("AUCell_buildRankings ...")
            # if (file.exists(file.path(dir_res, "aucell.AUCell_buildRankings.rds"))) {
            if (F) {
                cells_rankings <- read_rds(file.path(dir_res, "aucell.AUCell_buildRankings.rds"))
            } else {
                cells_rankings <- AUCell_buildRankings(
                    GetAssayData(xmo_sub, layer = "counts", assay = DefaultAssay(xmo_sub)),
                    plotStats = FALSE,
                    nCores = 10
                )
                write_rds(cells_rankings, file.path(dir_res, "aucell.AUCell_buildRankings.rds"))
            }
            cat("AUCell_calcAUC ...")
            # if (!file.exists(file.path(dir_res, 'aucell.AUCell_calcAUC.rds'))) {
            if (T) {
                cells_AUC <- AUCell_calcAUC(
                    geneSets = gene_set_use, rankings = cells_rankings,
                    aucMaxRank = ceiling(20 / 100 * nrow(cells_rankings)),
                    nCores = 10
                )
                ## 5% remove certain sets
                ## 20% over-claim gene sets
                ##
                write_rds(cells_AUC, file.path(dir_res, "aucell.AUCell_calcAUC.rds"))
            } else {
                cells_AUC <- read_rds(file.path(dir_res, "aucell.AUCell_calcAUC.rds"))
            }

            if (F) {
                ## not really needed
                cells_assignment <- AUCell_exploreThresholds(
                    cells_AUC,
                    plotHist = F, assignCells = T, nCores = 10
                )
                write_rds(cells_assignment, file.path(dir_res, "aucell.AUCell_exploreThresholds.rds"))
            }

            cat("getAUC ...")
            mat_auc <- t(getAUC(cells_AUC))
            mat_auc <- as.data.frame(mat_auc)
            if (!identical(Cells(xmo_sub), rownames(mat_auc))) {
                message("not matched cell names so fixe")
                mat_auc <- mat_auc[Cells(xmo_sub), ]
            }
            for (i in seq_along(names(gene_set_use))) {
                colnames(mat_auc)[i] <- name_module_set[i]
            }
            # mat_auc[1:3, 1:3]
            xmo_sub <- AddMetaData(xmo_sub, mat_auc)
        } else {
            #------ ERROR ------
            warning(c("unknown method: ", module_method))
        }

        df_module_sub <- xmo_sub@meta.data[, name_module_set]
        tmp_i <- match(Cells(xmo), rownames(df_module_sub))
        df_module <- df_module_sub[tmp_i, ]
        rownames(df_module) <- Cells(xmo)
        head(df_module_sub)
        head(df_module)
        rm(tmp_i)
        cat("Exporting... ")
        xmo <- AddMetaData(xmo, df_module)
        write_rds(
            xmo@meta.data[, name_module_set, drop = F],
            f_module
        )
        cat("[done]\n")
    }
}

#------------------ ~~~ Visualization ~~~ --------------------
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.visiumHD.R")
# source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")

theme_set(theme_pubr(base_size = 6, legend = "right"))
cli_h1("Visualization")
sp_xy_ratio <- get_visium_xy_ratio(xmo)
pdf_width <- 5
pdf_height <- 5
if (sp_xy_ratio > 1) {
    pdf_height <- pdf_height * sp_xy_ratio
} else {
    pdf_width <- pdf_width * sp_xy_ratio
}
print(c(pdf_width, pdf_height))

library(colorspace)
library(ggrastr)

#------ Wanted cells------
wanted_cells <- Cells(xmo)[xmo$celltype == cancer_ident]; str(wanted_cells)
p <- SpatialDimPlot(xmo, cells.highlight = wanted_cells, cols.highlight = c("red", NA)) +
    rremove("legend") + 
    labs(title = sprintf("Wanted %s cells (%s)", cancer_ident, comma(length(wanted_cells)))) +
    theme(aspect.ratio = sp_xy_ratio)
ggsave(
    filename = file.path(dir_res, "spatial.wanted_cells.pdf"),
    plot = p,
    width = pdf_width + 3,
    height = pdf_height, useDingbats = F
)
p <- DimPlot(
    xmo,
    cells.highlight = wanted_cells,
    cols.highlight = c("red", NA),
    label = F, alpha = 1
) +
    rremove("legend") +
    labs(title = sprintf("Wanted %s cells (%s)", cancer_ident, comma(length(wanted_cells)))) +
    theme(aspect.ratio = 1)
ggsave(
    filename = file.path(dir_res, "dimplot.wanted_cells.pdf"),
    plot = p, width = 4, height = 4, useDingbats = F
)
#------ Basic QC overview ------

cli_h2("Plotting nFeature and nCount")
xmo$log2nCount_Spatial <- log2(xmo$nCount_Spatial + 1)
for (y in c("nFeature_Spatial", "nCount_Spatial", "log2nCount_Spatial")) {
    p <- SpatialFeaturePlot(
        xmo,
        features = y,
        # pt.size.factor = 5,
        min.cutoff = "q1", max.cutoff = "q99",
        alpha = 1
    ) + theme(aspect.ratio = sp_xy_ratio)
    # p <- ggrastr::rasterise(p, layers='Spatial',dpi = 300) # points and the image are merged together
    ggsave(
        filename = file.path(dir_res, sprintf("spatial.feature.%s.pdf", y)),
        plot = p,
        width = pdf_width + 3,
        height = pdf_height, useDingbats = F
    )
    p <- FeaturePlot(
        xmo,
        features = y,
        min.cutoff = "q1", max.cutoff = "q99",
        label = T, alpha = 1, raster = T, pt.size = 3
    ) + theme(aspect.ratio = 1)
    ggsave(
        filename = file.path(dir_res, sprintf("dimplot.feature.%s.pdf", y)),
        plot = p,
        width = 4,
        height = 4, useDingbats = F
    )
}


# module_score_prefix <- paste0('module_score')
# module_score_prefix <- paste0('AMSREV_score')
# module_score_prefix <- paste0('ucell_score')
# module_score_prefix <- paste0('aucell_score')

for (module_method in method_opts) {
    cli_h2(module_method)
    module_score_prefix <- paste0(module_method, "_score")
    name_module_set <- paste0(module_score_prefix, "_", names(gene_set_use))

    vrange_module_set <- quantile(as.matrix(xmo@meta.data[, name_module_set]),
        c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = T
    )
    print(vrange_module_set)

    # img_pt_size <- 1.3
    # if (ncol(xmo) > 5000) {
    #     img_pt_size <- 1
    # }
    # if (ncol(xmo) > 10000) {
    #     img_pt_size <- 0.8
    # }
    # if (ncol(xmo) > 20000) {
    #     img_pt_size <- 0.5
    # }

    #------ viz module scores on spatial and umap ------
    cli_h3("viz module scores on spatial and umap")
    if (F) {
        ## this does not use the universal scale bar
        p1a <- lapply(name_module_set, function(mx) {
            SpatialFeaturePlot(
                xmo,
                features = mx,
                # pt.size.factor = 5,
                min.cutoff = vrange_module_set[1],
                max.cutoff = vrange_module_set[5],
                alpha = 1
            ) +
                theme(aspect.ratio = sp_xy_ratio) +
                labs(title = mx, color = mx) +
                scale_color_gradientn(
                    na.value = "#f8f8ff00",
                    colors = c(hcl.colors(n = 5, palette = "RdYlBu", rev = T)),
                    # values=rescale(vrange_module_set),
                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
                )
        })
        ggsave(
            file.path(dir_res, sprintf("%s.imgfeatureplot.birdview_gene_set.pdf", module_score_prefix)),
            ggarrange(
                plotlist = lapply(p1a, function(p) p + rremove("legend")),
                nrow = round(sqrt(length(name_module_set))), ncol = round(sqrt(length(name_module_set)))
            ),
            width = pdf_width * 1.5 * round(sqrt(length(name_module_set))),
            height = pdf_height * 1.5 * round(sqrt(length(name_module_set))), 
            useDingbats = F, limitsize = F
        )
        ggsave(
            file.path(dir_res, sprintf("%s.imgfeatureplot.birdview_gene_set.legend.pdf", module_score_prefix)),
            ggarrange(
                plotlist = lapply(p1a, function(x) ggpubr::get_legend(x) %>% as_ggplot()),
                nrow = round(sqrt(length(name_module_set))), ncol = round(sqrt(length(name_module_set)))
            ),
            width = 3, height = 3, 
            useDingbats = F, limitsize = F
        )
    }


    p1a <- SpatialFeaturePlot(xmo,
        features = name_module_set, keep.scale = 'all',
        ncol = round(sqrt(length(name_module_set)))
    ) & theme(aspect.ratio = sp_xy_ratio)
    ggsave(
        file.path(dir_res, sprintf("%s.imgfeatureplot.birdview_gene_set.pdf", module_score_prefix)),
        p1a,
        width = pdf_width * 1.5 * round(sqrt(length(name_module_set))),
        height = pdf_height * 1.5 * round(sqrt(length(name_module_set))) + 2, 
        useDingbats = F, limitsize = F
    )
    # ggsave(
    #         file.path(dir_res, sprintf("%s.imgfeatureplot.birdview_gene_set.legend.pdf", module_score_prefix)),
    #         ggpubr::get_legend(p1a) %>% as_ggplot()), 
    #         width = 3, height = 3, 
    #         useDingbats = F, limitsize = F
    #     )
    ggsave(
        file.path(dir_res, sprintf("%s.imgfeatureplot.birdview_gene_set.legend.pdf", module_score_prefix)),
        ggarrange(
            plotlist = lapply(p1a, function(x) ggpubr::get_legend(x) %>% as_ggplot()),
            nrow = round(sqrt(length(name_module_set))), ncol = round(sqrt(length(name_module_set)))
        ),
        width = 3, height = 3, 
        useDingbats = F, limitsize = F
    )
    
    p1b <- FeaturePlot(
        xmo,
        # cells = cancer_cnames,
        features = name_module_set,
        min.cutoff = vrange_module_set[1],
        max.cutoff = vrange_module_set[5],
        # keep.scale = 'all',
        # cols = c('lightgrey', 'maroon1'),
        pt.size = 3, ncol = round(sqrt(length(name_module_set))),
        raster = T, raster.dpi = c(1024, 1024)
    ) & theme(aspect.ratio = 1) & my_scatter_themevoid 
        # scale_color_gradientn(
        #     na.value = "ghostwhite",
        #     colors = c(hcl.colors(n = 5, palette = "RdYlBu", rev = T)),
        #     # values=rescale(vrange_module_set),
        #     guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
        # )

    # p1b
    ggsave(
        file.path(dir_res, sprintf("%s.featureplot.birdview_gene_set.pdf", module_score_prefix)),
        p1b ,
        width = 7.5 * round(sqrt(length(name_module_set))),
        height = 7.5 * round(sqrt(length(name_module_set)))
    )
    ggsave(
        file.path(dir_res, sprintf("%s.featureplot.birdview_gene_set.legend.pdf", module_score_prefix)),
        ggarrange(
            plotlist = lapply(p1b, function(x) ggpubr::get_legend(x) %>% as_ggplot()),
            nrow = round(sqrt(length(name_module_set))), ncol = round(sqrt(length(name_module_set)))
        ),
        width = 3, height = 3
    )


    #------ compare modules on cancer ident only ------
    cli_h3("compare modules on cancer ident only")

    tmp <- FetchData(xmo, cells = Cells(xmo)[Idents(xmo) %in% cancer_ident], vars = name_module_set)
    tmp <- rownames_to_column(tmp, "cellname")
    tmp <- tidyr::pivot_longer(tmp, !cellname, names_to = "module", values_to = "value")
    tmp$module_lab <- str_remove_all(tmp$module, pattern = module_score_prefix)
    all_pairwise_comparison_index <- function(n, m = 2) {
        choice_mat <- t(combn(x = 1:n, m = m))
        lapply(1:nrow(choice_mat), function(i) {
            choice_mat[i, ]
        })
    }
    tmp_name_module_lab <- sort(unique(tmp$module_lab))
    my_comparisons <- lapply(
        all_pairwise_comparison_index(length(tmp_name_module_lab)),
        function(i) tmp_name_module_lab[i]
    )
    p2 <- ggviolin(tmp, x = "module_lab", y = "value") +
        stat_mean() +
        scale_x_discrete(
            labels = ruok::pretty_table2str(formatC(tapply(tmp$value, tmp$module_lab, mean, na.rm = T), digits = 2))
        ) +
        stat_compare_means(comparisons = my_comparisons) +
        labs(y = module_score_prefix, title = sprintf("ident = %s", cancer_ident))

    ggsave(
        file.path(dir_res, sprintf("%s.vlnplot.birdview_over_cancer_ident.pdf", module_score_prefix)),
        p2,
        width = length(tmp_name_module_lab) * 1.2,
        height = 5, useDingbats = F, limitsize = F
    )

    ## Save the psbulk vector for any future use
    dir_res
    stats_report <- tmp %>% group_by(module) %>% 
        summarise(
            mean = mean(value, na.rm = T),
            median = median(value, na.rm = T),
            sd = sd(value, na.rm = T), 
            SDE = stderror(value)
        )
    write_csv(stats_report, file.path(dir_res, sprintf("stats_report.%s.%s.csv", gene_suit_name, module_score_prefix)))
    stats_report[, c(1,2), drop=F] %>% deframe() %>% 
        write_rds(file.path(dir_res, sprintf("psbulk_vector.%s.%s.rds", gene_suit_name, module_score_prefix)))
    cat("[done]\n")
}

cli_alert_success("DONE")
