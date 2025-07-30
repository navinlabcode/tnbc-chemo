suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Add module scores
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Yun Yan (yun.yan@uth.tmc.edu)
# Date Created: 2024-08-17
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
    library(fs)
    library(Signac)
    library(clustree)
    library(UCell)
    library(AUCell)
    library(BiocParallel)
    my_scatter_themevoid <- theme_pubr(base_size = 6, legend = "right") %+replace% theme(
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = rel(1)),
        axis.line = element_blank()
    )
    source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
    source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/util.nmf.viz.R")
})
options <- commandArgs(trailingOnly = TRUE)

# For example: HBCA lineage for singel-cell data
if (F) {
    f_in <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas",
        "rds_rna-integrate/pat102/lv01.aneuploidy_tri_type.aneuploid.pure5",
        "ready.sr3.rds"
    )
    f_gene_set <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/meta/hbca_epithelial_lineage_siyuan/epithelial_lineage.rds"
    gene_suit_name <- "hbca_epithelial_lineage_siyuan_head15"
    f_cellmeta <- file.path(
        "/volumes/USR1/yyan/project/tnbc_pre_atlas",
        "rds_rna-integrate/pat102/lv01.aneuploidy_tri_type.aneuploid.pure5",
        "sr3_metadata.df.rds"
    )
    assay_use <- "RNA"
    ident_by_what <- "patient"
    use_head_n_genes <- 15
}

if (length(options) > 0) {
    f_in <- options[1]
    f_gene_set <- options[2]
    gene_suit_name <- options[3]
    assay_use <- options[4]
    f_cellmeta <- options[5]
    ident_by_what <- options[6]
}

#------------------ ~~~ Readin ~~~ --------------------
cli_h1("Readin")
message(f_in)
sr3 <- read_rds(f_in)
DefaultAssay(sr3) <- assay_use
print(sr3)

table(sr3@meta.data[[ident_by_what]], useNA = "always") %>% print()
method_opts <- c("module", "AMSREV", "ucell", "aucell")
method_opts <- c("module", "ucell", "aucell")
method_opts <- c("module", "aucell")
# method_opts <- c("module")
cli_li(method_opts)
#------------------ ~~~ gene lists ~~~ --------------------
gene_set_use <- readRDS(f_gene_set)
gene_set_use <- lapply(gene_set_use, function(v) {
    intersect(v, rownames(sr3@assays[[DefaultAssay(sr3)]]@counts))
})

# gene_set_use <- lapply(gene_set_use, rand_slice, n = 50)
gene_set_use <- lapply(gene_set_use, head, n= use_head_n_genes)
str(gene_set_use)
#------------------ ~~~ setup output dir ~~~ --------------------
dir_res <- file.path(dirname(f_in), sprintf("modulescore_%s", gene_suit_name))
dir_create(dir_res)
write_rds(gene_set_use, file.path(dir_res, "gene_set.list.rds"))
write_gmx2(gene_set_use, file.path(dir_res, "gene_set.gmx.csv"))
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
        cat("Load the existing results... ")

        if (any(str_detect(colnames(sr3@meta.data), module_score_prefix))) {
            cat("the object had results so overwrite...")
        }

        df_module <- read_rds(f_module)
        stopifnot(identical(rownames(df_module), Cells(sr3)))
        sr3 <- AddMetaData(sr3, df_module)
    } else {
        cat("Computing...")
        run_module_func <- switch(module_method,
            module = AddModuleScore,
            AMSREV = AMS_rev,
            ucell = AddModuleScore_UCell
        )
        module_score_prefix <- paste0(module_method, "_score")
        name_module_set <- paste0(module_score_prefix, "_", names(gene_set_use))

        if (module_method %in% c("module", "AMSREV")) {
            sr3 <- run_module_func(sr3, gene_set_use,
                name = module_score_prefix
            )
            for (i in seq_along(names(gene_set_use))) {
                i_old <- paste0(module_score_prefix, i)
                i_new <- name_module_set[i]
                sr3@meta.data[[i_new]] <- sr3@meta.data[[i_old]]
                sr3@meta.data[[i_old]] <- NULL
            }
        } else if (module_method == "ucell") {
            sr3 <- run_module_func(sr3, gene_set_use, name = "_UCell")
            for (i in seq_along(names(gene_set_use))) {
                i_old <- paste0(names(gene_set_use)[i], "_UCell")
                i_new <- name_module_set[i]
                sr3@meta.data[[i_new]] <- sr3@meta.data[[i_old]]
                sr3@meta.data[[i_old]] <- NULL
            }
        } else if (module_method == "aucell") {
            cat("AUCell_buildRankings ...")
            if (file.exists(file.path(dir_res, "aucell.AUCell_buildRankings.rds"))) {
                cells_rankings <- read_rds(file.path(dir_res, "aucell.AUCell_buildRankings.rds"))
            } else {
                cells_rankings <- AUCell_buildRankings(
                    GetAssayData(sr3, slot = "counts", assay = DefaultAssay(sr3)),
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
            if (!identical(Cells(sr3), rownames(mat_auc))) {
                message("not matched cell names so fixe")
                mat_auc <- mat_auc[Cells(sr3), ]
            }
            for (i in seq_along(names(gene_set_use))) {
                colnames(mat_auc)[i] <- name_module_set[i]
            }
            # mat_auc[1:3, 1:3]
            sr3 <- AddMetaData(sr3, mat_auc)
        } else {
            warning(c("unknown method: ", module_method))
        }

        cat("Exporting... ")
        write_rds(sr3@meta.data[, name_module_set, drop = F], f_module)
        cat("[done]\n")
    }
}

#------------------ ~~~ Visualization on UMAP ~~~ --------------------

theme_set(theme_pubr(base_size = 7, legend = "right"))
cli_h1("Visualization")

library(colorspace)

if ("umap" %in% names(sr3@reductions)) {
    for (feature_i in paste0(c("nCount_", "nFeature_"), assay_use)) {
        message(feature_i)
        p <- FeaturePlot(
            sr3,
            features = feature_i, pt.size = 3,
            reduction = "umap",
            max.cutoff = "q95",
            min.cutoff = "q5", order = T, raster = T, raster.dpi = c(1024, 1024),
            cols = c("lightgrey", "maroon1")
        ) +
            my_scatter_themevoid
        ggsave(file.path(dir_res, sprintf("featureplot.%s.pdf", feature_i)),
            p,
            width = 7, height = 7, useDingbats = F
        )
    }

    for (module_method in method_opts) {
        cli_h2(module_method)
        module_score_prefix <- paste0(module_method, "_score")
        name_module_set <- paste0(module_score_prefix, "_", names(gene_set_use))
        length(name_module_set)
        vrange_module_set <- quantile(
            as.matrix(sr3@meta.data[, name_module_set]),
            c(0.05, 0.25, 0.5, 0.75, 0.95)
        )
        print(vrange_module_set)

        p1b <- FeaturePlot(
            sr3,
            features = name_module_set,
            min.cutoff = vrange_module_set[1],
            max.cutoff = vrange_module_set[5],
            # cols = c('lightgrey', 'maroon1'),
            pt.size = 3,
            raster = T, raster.dpi = c(1024, 1024)
        ) &
            # scale_color_viridis_c(option = 'H', guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)) &
            scale_color_gradientn(
                colors = c(hcl.colors(n = 5, palette = "RdYlBu", rev = T)),
                # values=rescale(vrange_module_set),
                guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
            ) &
            labs(color = "") & theme(aspect.ratio = 1) & my_scatter_themevoid

        # p1b
        ggsave(
            file.path(dir_res, sprintf("%s.featureplot.birdview_gene_set.pdf", module_score_prefix)),
            p1b & rremove("legend"),
            width = sqrt(length(name_module_set)) * 1.875,
            height = sqrt(length(name_module_set)) * 1.875,
            useDingbats = F
        )
        ggsave(
            file.path(dir_res, sprintf("%s.featureplot.birdview_gene_set.legend.pdf", module_score_prefix)),
            as_ggplot(ggpubr::get_legend(p1b)),
            width = 1.5, height = 1.5, useDingbats = F
        )
    }
}
#------------------ ~~~ Average heatmap of modules ~~~ --------------------
if (is.null(f_cellmeta)) {
    message("No cell metadata file provided so skip Average heatmap")
    q()
}
if (is.null(ident_by_what)) {
    message("No ident_by_what provided so skip Average heatmap")
    q()
}

df_cellmeta <- read_rds(f_cellmeta)
cli_alert_info(ident_by_what)
stopifnot(identical(rownames(df_cellmeta), Cells(sr3)))

if (!ident_by_what %in% colnames(sr3@meta.data)) {
    message("Add cell metadata to Seurat object")
    sr3 <- AddMetaData(sr3, df_cellmeta[, ident_by_what, drop = F])
}

cli_h1("Average heatmap of modules")

for (module_method in method_opts) {
    cli_h2(module_method)
    module_score_prefix <- paste0(module_method, "_score")
    name_module_set <- paste0(module_score_prefix, "_", names(gene_set_use))
    length(name_module_set)

    df_val <- sr3@meta.data[, c(ident_by_what, name_module_set), drop = F]
    head(df_val)

    df_val <- df_val %>%
        group_by(!!sym(ident_by_what)) %>%
        summarise(across(all_of(name_module_set), mean)) %>%
        ungroup() %>%
        column_to_rownames(ident_by_what)
    colnames(df_val) <- gsub(paste0(module_score_prefix, "_"), "", colnames(df_val))
    all(names(gene_set_use) %in% colnames(df_val))
    df_val <- df_val[, names(gene_set_use)]
    mat_val <- as.matrix(df_val) # idents x modules
    write_rds(mat_val, file.path(dir_res, sprintf("%s.avg_%s.orig.rds", module_score_prefix, ident_by_what)))
    dim(mat_val)
    for (val_type in c("orig", "zscore", "rescale")) {
        if (val_type == "zscore") {
            mat_viz <- scale(mat_val)
        } else if (val_type == "orig") {
            mat_viz <- mat_val
        } else if (val_type == "rescale") {
            mat_viz <- scale_nmf_gene_loading(mat_val)
        } else {
            stop("unknown val_type")
        }
        stopifnot(identical(dim(mat_viz), dim(mat_val)))
        mat_viz <- t(mat_viz)
        vrange <- quantile(mat_viz, c(0.05, 0.95))
        if (val_type %in% c("orig", "zscore")) {
            col_use <- heatmap_color_fun_cont(vrange[1], vrange[2], "Blue-Red 2")
        } else {
            col_use <- heatmap_color_fun_cont(vrange[1], vrange[2], "Viridis")
        }
        col_use <- heatmap_color_fun_cont(vrange[1], vrange[2], "Blue-Red 2")
        Heatmap(
            mat_viz,
            col = col_use,
            show_row_names = T, row_names_side = "right",
            show_column_names = T, column_names_side = "top",
            cluster_rows = F,
            cluster_columns = T,
            column_title = "Module",
            row_title = ident_by_what,
            name = val_type,
            column_dend_height = unit(0.1, "inch"),
            height = unit(5, "inch"),
            width = unit(5 * (ncol(mat_viz) / nrow(mat_viz)), "inch"),
        ) -> p
        pdf(file.path(dir_res, sprintf("heatmap.%s.%s.avg_by_%s.pdf", module_score_prefix, val_type, ident_by_what)),
            height = 5 + 3,
            width = 5 * (ncol(mat_viz) / nrow(mat_viz)) + 5
        )
        draw(p)
        dev.off()
    }
}
cli_alert_success("DONE")
