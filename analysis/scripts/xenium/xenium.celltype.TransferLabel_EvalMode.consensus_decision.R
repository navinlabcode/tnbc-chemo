suppressPackageStartupMessages(library(magrittr))
"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal: Consensus determine cell identity (`cat_transfer`)
# using the 10 portion genes results.
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
  library(ruok)
  library(fs)
  library(stringr)
  library(Seurat) 
  library(pbapply)
  source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.xenium.R")
  source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.pal.R")
  source("~/project/tumor_plasticity/sandbox/pre_atlas_codes/util.stats.R")
  library(nanoparquet)
  library(ComplexHeatmap)
  library(glue)
})
library(future)
options(future.globals.maxSize = 10 * 1024^3)
cmdargs <- commandArgs(trailingOnly = TRUE)
print(cmdargs)


if (length(cmdargs) > 0) {
  f_sp <- cmdargs[[1]]
  cat_transfer <- cmdargs[[2]] # celltype | cell_state_paper
  # sc_assay           <- cmdargs[[3]] # RNA
  # sp_assay           <- cmdargs[[4]] # Xenium
} else {
  f_sp <- "/volumes/USR1/yyan/project/tnbc_xenium/data/ART94/nonbinarized_pca/xenium_nonbinarized_pca.seurat.rds"
  sp_assay <- "Xenium"
  f_sc <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/downto_1000_by_cell_state_paper/ready.sr3.rds"
  cat_transfer <- "celltype" ## can change to 'cell_state_paper'.
}
if (T) {
  sc_assay <- "RNA"
  sp_assay <- "Xenium"
}

#------------------ ~~~ Fixed parameters ~~~ --------------------
n_gene_portions <- 10 ## do not change
str_gene_portion_folder_path_prefix <- "TransferLabelEval.ref_RNA.query_Xenium.gene_portion_" ## do not change

#------------------ ~~~ Readin ~~~ --------------------
# sc <- read_rds(f_sc); DefaultAssay(sc) <- sc_assay
sp <- read_rds(f_sp)
DefaultAssay(sp) <- sp_assay




#------------------ ~~~ Load results of the 10-gene-portions ~~~ --------------------

dir_portions_list <- file.path(
  dirname(f_sp),
  sprintf(
    "TransferLabelEval.ref_%s.query_%s.gene_portion_%s",
    sc_assay, sp_assay, 1:n_gene_portions
  )
)

dir_portion_transfer_list <- file.path(dir_portions_list, cat_transfer)
f_transfer_list <- file.path(dir_portion_transfer_list, "TransferData_celltype.dataframe.rds")
names(f_transfer_list) <- as.character(1:n_gene_portions)
all(file.exists(f_transfer_list))

f_transfer_list <- f_transfer_list[file.exists(f_transfer_list)]
cli_li(f_transfer_list)

df_list <- lapply(f_transfer_list, function(x) {
  o <- read_rds(x)
  o <- rownames_to_column(o, "cellname")
  o
})

print(lapply(df_list, head, n = 5))

#------------------ ~~~ Set up output dir ~~~ --------------------

dir_res <- file.path(
  dirname(f_sp),
  sprintf(
    "TransferLabelEval.ref_%s.query_%s.Consensus",
    sc_assay, sp_assay
  ),
  cat_transfer
)
fs::dir_create(dir_res)

#------------------ ~~~ Setup to help visualization ~~~ --------------------
shared_column_names <- Reduce(union, sapply(df_list, colnames))
cat_levels <- shared_column_names[str_detect(shared_column_names, "prediction.score.")] %>%
  str_remove_all(., pattern = "prediction.score.")

if (cat_transfer %in% c("celltype")) {
  pal_cat <- pal_celltypes
}
cat_levels <- intersect(names(pal_cat), cat_levels)
print(cat_levels)

if (T) {
  pl <- enframe(pal_cat) %>%
    mutate(name = factor(name, levels = name)) %>%
    ggplot() +
    geom_point(aes(x = name, y = factor(1), color = name), size = 3) +
    scale_color_manual(values = pal_cat) +
    labs(color = cat_transfer)
  ggsave(file.path(dir_res, sprintf("legend.%s.pdf", cat_transfer)),
    as_ggplot(get_legend(pl)),
    width = 3, height = 3, useDingbats = F
  )
  rm(pl)
}

#------------------ ~~~ Whether results are stable across protions ~~~ --------------------
# a long df: cell | portion_id | predicted.id
df_cell_portion_ident <- lapply(names(df_list), function(x) {
  df <- df_list[[x]]
  df$portion <- x
  return(df[, c("cellname", "portion", "predicted.id")])
}) %>% do.call("rbind", .)
head(df_cell_portion_ident, 3)
#     cellname portion predicted.id
# 1 aaaccmkb-1       1        Fibro
df_cell_portion_ident$predicted.id <- factor(
  as.character(df_cell_portion_ident$predicted.id),
  levels = cat_levels
)
df_cell_portion_ident$portion <- as.factor(as.numeric(df_cell_portion_ident$portion))
p <- df_cell_portion_ident %>%
  ggplot(., aes(x = portion, fill = predicted.id)) +
  geom_bar() +
  scale_fill_manual(values = pal_cat) +
  labs(y = "ncells")
p <- p + labs(caption = report_chisq(
  chisq.test(table(
    df_cell_portion_ident$predicted.id,
    df_cell_portion_ident$portion
  ))
))
# p
ruok::ggsave2(file.path(dir_res, "barplot.ident_by_portions"),
  p + rremove("legend"),
  width = 4, height = 3
)

#------------------ ~~~ Determine the consensus idents ~~~ --------------------
## a wide data.frame: cell | css_id | nportion_css_id
df_css_ident <- df_cell_portion_ident %>%
  dplyr::group_by(cellname, predicted.id) %>%
  dplyr::summarise(nportion = n())
df_css_ident <- df_css_ident %>%
  dplyr::group_by(cellname) %>%
  dplyr::slice(nnet::which.is.max(nportion))
colnames(df_css_ident) <- c("cellname", "css_id", "nportion.css_id")
print(head(df_css_ident)) # ! export
#   cellname   css_id nportion.css_id
# 1 aaaccmkb-1 Fibro                5
print(table(df_css_ident$css_id))

#------------------ ~~~ Prepare diagnosis of the consensus idents ~~~ --------------------
## a wide data.frame: cell | nportion.<EACH_CAT>
head(df_list[[1]], 3)
dflong_ident_nportion <- df_cell_portion_ident %>%
  dplyr::group_by(cellname, predicted.id) %>%
  dplyr::summarise(nportion = n())
head(dflong_ident_nportion)
df_cell_nportioncat <- dflong_ident_nportion %>%
  pivot_wider(
    id_cols = "cellname",
    names_from = "predicted.id",
    names_prefix = "nportion.",
    values_from = "nportion",
    values_fill = 0
  )
head(df_cell_nportioncat, 1) # ! export

## a wide data.frame: cell | prediction.score.<EACH_CAT>
df_cell_portion_probcat <- lapply(names(df_list), function(x) {
  df <- df_list[[x]]
  df$portion <- x
  col_str_use <- paste0("prediction.score.", cat_levels)
  if (!all(col_str_use %in% colnames(df))) {
    missing_col_str <- setdiff(col_str_use, colnames(df))
    for (tmp in missing_col_str) {
      df[[tmp]] <- 0
    }
    rm(tmp)
  }

  return(df[, c("cellname", "portion", col_str_use)])
}) %>% do.call("rbind", .)
head(df_cell_portion_probcat, 1)
#     cellname portion prediction.score.Tumor prediction.score.Mye
# 1 aaaccmkb-1       1             0.05885438            0.3650900
#   prediction.score.T prediction.score.B prediction.score.Fibro
# 1         0.06211791         0.04690388             0.42466257
#   prediction.score.Endo prediction.score.Peri
# 1            0.01571327            0.02665803
df_cell_meanprobcat <- df_cell_portion_probcat %>%
  dplyr::select(-portion) %>%
  dplyr::group_by(cellname) %>%
  dplyr::summarise_all(., mean, na.rm = T)
head(df_cell_meanprobcat) ## ! export

stopifnot(identical(sort(df_css_ident$cellname), sort(df_cell_meanprobcat$cellname)))
stopifnot(identical(sort(df_css_ident$cellname), sort(df_cell_nportioncat$cellname)))


df_css_res <- dplyr::left_join(df_css_ident, df_cell_nportioncat)
df_css_res <- dplyr::left_join(df_css_res, df_cell_meanprobcat)
# view(head(df_css_res))

df_css_res$prediction.score.css_id <- pbsapply(
  1:nrow(df_css_res), function(i) {
    x <- df_css_res$css_id[i]
    tmp <- paste0("prediction.score.", x)
    as.numeric(df_css_res[i, tmp])
  }
)

head(df_css_res)

# stopifnot(all(df_css_res$cellname %in% Cells(sp)))
if (!identical(df_css_res$cellname, Cells(sp))) {
  cli_alert_warning("reorder cell names to match with seurat object")
  tmp <- match(Cells(sp), df_css_res)
  df_css_res <- df_css_res[tmp, ]
  df_css_res$cellname <- Cells(sp)
}


if (T) {
  ## export
  write_rds(
    df_css_res,
    file.path(dir_res, sprintf("Consensus_TransferData_%s.dataframe.rds", "css_id"))
  )
  write_csv(
    df_css_res,
    file.path(dir_res, sprintf("Consensus_TransferData_%s.dataframe.csv", "css_id"))
  )
  nanoparquet::write_parquet(
    df_css_res,
    file.path(dir_res, sprintf("Consensus_TransferData_%s.dataframe.parquet", "css_id"))
  )
}
#------------------ ~~~ Diagnosis consensus idents ~~~ --------------------
#
# df_css_res <- read_rds(
#   file.path(dir_res, sprintf('Consensus_TransferData_%s.dataframe.rds', 'css_id')))
head(df_css_res)
df_css_res <- df_css_res %>% column_to_rownames("cellname")
# 1/3 heatmap: ask if a cell is specific to css_id. Comparing 1 vs the rest.
adhoc_heatmap <- function(
    df, ident_col,
    ident_levels = NULL,
    measure_class = c("nportion", "prediction.score"),
    pal_idents = NULL,
    n_downsample = 100, seed = 1026,
    ...) {
  measure_class <- match.arg(measure_class)
  ident_val_col <- paste0(measure_class, ".", ident_col)

  if (is.null(ident_levels)) {
    ident_levels <- unique(as.character(df[[ident_col]]))
  }
  if (is.null(pal_idents)) {
    pal_idents <- structure(rainbow(n = length(ident_levels)), names = ident_levels)
  }

  value_cols <- paste0(measure_class, ".", ident_levels)

  set.seed(seed)
  df <- df %>%
    dplyr::group_by_at(ident_col) %>%
    dplyr::sample_n(size = n_downsample)
  df[["tmp"]] <- df[[ident_val_col]] * -1
  df <- df %>% dplyr::arrange(tmp, .by_group = T)

  mat <- df[, value_cols] %>% as.matrix()
  colnames(mat) <- str_remove_all(colnames(mat), measure_class)

  ComplexHeatmap::Heatmap(
    matrix = mat, name = measure_class,
    cluster_rows = F, cluster_columns = F,
    column_names_side = "top",
    row_title = sprintf("random %d cells per identity", n_downsample),
    left_annotation = rowAnnotation(
      cluster = df[[ident_col]],
      col = list(cluster = pal_idents)
    ),
    use_raster = T, raster_by_magick = T,
    ...
  )
}


pdf(
  file.path(
    dir_res,
    sprintf("diagnose.css_id.%s.1vsothers.heatmap.pdf", "nportion")
  ),
  width = 3.5, height = 5, useDingbats = F
)
p <- adhoc_heatmap(df_css_res,
  ident_col = "css_id",
  ident_levels = cat_levels,
  measure_class = "nportion",
  pal_idents = pal_cat
)
draw(p)
dev.off()
pdf(
  file.path(
    dir_res,
    sprintf("diagnose.css_id.%s.1vsothers.heatmap.pdf", "prediction.score")
  ),
  width = 3.5, height = 5, useDingbats = F
)
p <- adhoc_heatmap(df_css_res,
  ident_col = "css_id",
  ident_levels = cat_levels,
  measure_class = "prediction.score",
  pal_idents = pal_cat
)
draw(p)
dev.off()

# barplot: focus on each css_id and check if all 10 portions get the same result
#------ nportion split by nportion.css_id ------
p1 <- table(df_css_res$nportion.css_id) %>%
  as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) +
  geom_col() +
  geom_text(aes(y = Freq, label = comma(Freq)), vjust = -.5, size = 6 / .pt) +
  scale_x_discrete(limits = factor(1:10)) +
  labs(y = "num cells", x = "num portions having the css_id")
# p1
p2 <- prop.table(table(df_css_res$nportion.css_id)) %>%
  as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) +
  geom_col() +
  geom_text(aes(y = Freq, label = percent(signif(Freq, digits = 2))), vjust = -.5, size = 6 / .pt) +
  scale_x_discrete(limits = factor(1:10)) +
  labs(y = "% cells", x = "num portions having the css_id") +
  scale_y_continuous(labels = percent)
# p2
p <- (p1 / p2)
ggsave2(file.path(dir_res, "barplot.check_sensitivity.css_id.nportion"), p, width = 4, height = 4)

pdf(file.path(dir_res, "barplot.check_sensitivity.splitby_css_id.nportion.pdf"),
  width = 4, height = 4.3, onefile = T, useDingbats = F
)
for (cat_lv in cat_levels) {
  message(cat_lv)
  p1 <- table(df_css_res[df_css_res$css_id == cat_lv, "nportion.css_id"]) %>%
    as.data.frame() %>%
    ggplot(aes(x = Var1, y = Freq)) +
    geom_col(fill = pal_cat[cat_lv]) +
    geom_text(aes(y = Freq, label = comma(Freq)), vjust = -.5, size = 6 / .pt) +
    scale_x_discrete(limits = factor(1:10)) +
    labs(y = sprintf("num of %s cells", cat_lv), x = "num portions having the css_id")
  # p1
  p2 <- table(df_css_res[df_css_res$css_id == cat_lv, "nportion.css_id"]) %>%
    prop.table() %>%
    as.data.frame() %>%
    ggplot(aes(x = Var1, y = Freq)) +
    geom_col(fill = pal_cat[cat_lv]) +
    geom_text(aes(y = Freq, label = percent(signif(Freq, digits = 2))), vjust = -.5, size = 6 / .pt) +
    labs(y = sprintf("%% of %s cells", cat_lv), x = "num portions having the css_id") +
    scale_y_continuous(labels = percent) +
    scale_x_discrete(limits = factor(1:10))
  # p2

  p <- (p1 / p2)
  print(p)
}
dev.off()

#------ prob split by nportion.css_id ------
p <- df_css_res %>%
  ggplot(., aes(x = as.factor(nportion.css_id), y = prediction.score.css_id)) +
  geom_violin(scale = "area", fill = "lightgrey", color = NA) +
  stat_mean(color = "black", pch = 16) +
  scale_x_discrete(limits = factor(1:10)) +
  labs(y = sprintf("mean prob"), x = "num portions having the css_id")

ggsave2(file.path(dir_res, "violin.check_sensitivity.css_id.prob"), p, width = 4, height = 2.5)

pdf(file.path(dir_res, "violin.check_sensitivity.splitby_css_id.prob.pdf"),
  width = 4, height = 2.5, onefile = T, useDingbats = F
)
for (cat_lv in cat_levels) {
  message(cat_lv)
  p <- df_css_res %>%
    dplyr::filter(css_id == cat_lv) %>%
    ggplot(., aes(x = as.factor(nportion.css_id), y = prediction.score.css_id)) +
    geom_violin(scale = "area", fill = pal_cat[cat_lv], color = NA) +
    stat_mean(color = "black", pch = 16) +
    scale_x_discrete(limits = factor(1:10)) +
    labs(y = sprintf("mean prob of being %s", cat_lv), x = "num portions having the css_id")
  print(p)
}
dev.off()


#------------------ ~~~ nGenes/nUMIs split by nportion.css_id ~~~ --------------------
## similar analysis as the prob split by nportion.css_id
df_css_res$cellname <- rownames(df_css_res)
df_css_res$nCount_Xenium <- sp@meta.data[df_css_res$cellname, "nCount_Xenium"]
df_css_res$nFeature_Xenium <- sp@meta.data[df_css_res$cellname, "nFeature_Xenium"]

for (basic_metric in c("nCount_Xenium", "nFeature_Xenium")) {
  # print(mean(df_css_res[[basic_metric]]))
  p <- df_css_res %>%
    dplyr::mutate(nportion.css_id = as.factor(nportion.css_id)) %>%
    ggplot(., aes_string(x = "nportion.css_id", y = basic_metric)) +
    geom_violin(scale = "area", fill = "lightgrey", color = NA) +
    stat_mean(color = "black", pch = 16) +
    scale_x_discrete(limits = factor(1:10)) +
    labs(y = sprintf(basic_metric), x = "num portions having the css_id")

  p <- p + scale_y_log10() + annotation_logticks(sides = "l")

  if (basic_metric == "nFeature_Xenium") {
    p <- p + geom_hline(yintercept = 500, lty = "dashed")
  }
  ggsave2(file.path(dir_res, glue("violin.check_sensitivity.css_id.{basic_metric}")),
    p,
    width = 4, height = 2.5
  )

  pdf(
    file.path(
      dir_res,
      glue("violin.check_sensitivity.splitby_css_id.{basic_metric}.pdf")
    ),
    width = 4, height = 2.5, onefile = T, useDingbats = F
  )
  for (cat_lv in cat_levels) {
    message(cat_lv)
    p <- df_css_res %>%
      dplyr::filter(css_id == cat_lv) %>%
      dplyr::mutate(nportion.css_id = as.factor(nportion.css_id)) %>%
      ggplot(., aes_string(x = "nportion.css_id", y = basic_metric)) +
      geom_violin(scale = "area", fill = pal_cat[cat_lv], color = NA) +
      stat_mean(color = "black", pch = 16) +
      scale_x_discrete(limits = factor(1:10)) +
      labs(
        y = sprintf("%s of being %s", basic_metric, cat_lv),
        x = "num portions having the css_id"
      )
    p <- p + scale_y_log10() + annotation_logticks(sides = "l")
    if (basic_metric == "nFeature_Xenium") {
      p <- p + geom_hline(yintercept = 500, lty = "dashed")
    }
    print(p)
  }
  dev.off()
}

# try(df_css_res$nFeature_Xenium <- NULL)
# try(df_css_res$nCount_Xenium <- NULL)


#------------------ ~~~ Propose the final identity ~~~ --------------------
css_id_full <- as.character(df_css_res$css_id)
css_id_full_prob <- df_css_res$prediction.score.css_id
is_bad_cell <- df_css_res$nportion.css_id != 10
table(is_bad_cell)
# is_bad_cell <- df_css_res$nportion.css_id < 9; table(is_bad_cell)

css_id_full[is_bad_cell] <- "LOWCONF"
css_id_full <- factor(css_id_full, levels = c(cat_levels, "LOWCONF"))
table(css_id_full)
df_css_res$css_id_full <- css_id_full
df_css_res$prediction.score.css_id_full <- css_id_full_prob

df_css_res$cellname <- rownames(df_css_res)
head(df_css_res)


if (T) {
  ## export
  write_rds(
    df_css_res,
    file.path(dir_res, sprintf("Consensus_TransferData_%s.dataframe.rds", "css_id_full"))
  )
  write_csv(
    df_css_res,
    file.path(dir_res, sprintf("Consensus_TransferData_%s.dataframe.csv", "css_id_full"))
  )
  nanoparquet::write_parquet(
    df_css_res,
    file.path(dir_res, sprintf("Consensus_TransferData_%s.dataframe.parquet", "css_id_full"))
  )
}


#------------------ ~~~ Export transferred label to Xenium Explorer ~~~ --------------------
for (ident_str in c("css_id", "css_id_full")) {
  df_css_res %>%
    df_to_exnium_explorer(., "cellname", ident_str) %>%
    write_csv(., file.path(
      dir_res,
      sprintf("to_xenium_explorer_groups.%s.csv", ident_str)
    ))
}

#------------------ ~~~ Spatial visualizing the consensus idents ~~~ --------------------
# use snippet code
stopifnot(identical(rownames(df_css_res), Cells(sp)))
dir_snippet_viz <- dir_res
xmo <- AddMetaData(sp, df_css_res)
print(xmo)
pal_z <- pal_cat
for (viz_what in c("css_id", "css_id_full")) {
  source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_categorical.R")
}
table(Idents(xmo))
xmo$css_id_clean <- xmo$css_id_full
xmo <- subset(xmo, cells = Cells(xmo)[xmo$css_id_full != "LOWCONF"])
viz_what <- "css_id_clean"
source("/volumes/USR1/yyan/project/tumor_plasticity/sandbox/pre_atlas_codes/xenium.snippet_viz_categorical.R")


cat("[done] R script")
timestamp()

#------------------ ~~~ Visualize - transferred genes ~~~ --------------------
# to-do in a separate script
# correlation scatter plot of psbulk. Each dot is a gene.
