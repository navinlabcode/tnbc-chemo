#---------------------------
# Response-related in ARTEMIS              ----
#---------------------------

#------ Initiate the data ------

# df_surv: data frame pCR_status, archetypes
# mat: bulk data
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(pROC)
if (T) {
  # psbulk.prepare.R
  dir_lib <- "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/psbulk/"
  df_surv <- read_rds(file.path(dir_lib, "patient_metainfo.rds"))
  mat <- read_rds(file.path(dir_lib, "normalized.matrix.rds"))
  print(range(mat, na.rm = T)) # 0 18.3432
} else {
  # use psbulk.prepare.R
}

rownames(df_surv) <- df_surv$patient
pat_use <- intersect(colnames(mat), df_surv$patient)
str(pat_use)
colnames(mat)
mat <- mat[, pat_use]
df_surv <- df_surv[pat_use, ]
stopifnot(all.equal(df_surv$patient, colnames(mat)))

dir_proj <- file.path(
  "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102",
  "lv01.aneuploidy_tri_type.aneuploid.pure5", "survival_analysis",
  "artemis_chemoyes.logistic_informative_genes"
)
fs::dir_create(dir_proj)


#------ gene list ------
dict_update_cellstatenames <- read_tsv("/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/color.css", col_names = c("old", "color", "new"), col_types = c("c", "c", "c"))
dict_update_cellstatenames <- dict_update_cellstatenames %>% filter(!new %in% c(NA))
dict_update_cellstatenames <- dict_update_cellstatenames[, c("old", "new")] %>% deframe()

if (T) {
  module_content <- read_rds(
    file.path(
      "/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102",
      "lv01.aneuploidy_tri_type.aneuploid.pure5",
      "metamodule_fnmf",
      "MM_alt_clean_byscore_wardD2", "deliver.mm_markers.rds"
    )
  )
  str(module_content)
  module_focuse <- setdiff(names(module_content), "M14")
  module_focuse_bio <- c(
    "G2/M", "Mitochondria", "Ribosome", "Stress", "Interferon",
    "HLA", "S/G1", "Hypoxia", "Basal", "EpithelialDiff",
    "LumSec", "Cholestero", "StressER"
  )
  names(module_focuse_bio) <- module_focuse
  ## adhoc choice!!!
  publish_top_markers <- 30 # use all
  module_content_use <- lapply(module_content, function(x) head(x, publish_top_markers))
  module_content_use <- module_content_use[module_focuse]
  str(module_content_use)
  names(module_content_use) <- paste0(
    names(module_content_use), "__",
    make.names(module_focuse_bio[names(module_content_use)])
  )
  # db <- c(module_content_use, db_msigdb_H)
  str(module_content_use)

  nmf_genes <- read_rds("~/project/tnbc_pre_atlas//rds_rna-integrate/pat102/lv01.aneuploidy_tri_type.aneuploid.pure5/psbulk/fastnmf/rank4/deliver.nmf_markers.top.rds")
  # nmf_genes <- read_rds('~/project/tnbc_pre_atlas//rds_rna-integrate/pat102/lv01.aneuploidy_tri_type.aneuploid.pure5/psbulk/fastnmf/rank4/deliver.nmf_markers.rds')
  nmf_genes <- lapply(nmf_genes, head, n = 100)
  # str(nmf_genes)
  names(nmf_genes) <- str_replace(names(nmf_genes), "fNMF", "ARC")
  str(nmf_genes)
  ## cell state specific pCR/RD genes
  cellstate_pcr_rd_deg <- read_rds("~/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/findallmarker.ecotrait_feature.RD_vs_pCR/deg_collect.rds")
  cellstate_pcr_rd_deg$cluster <- dict_update_cellstatenames[cellstate_pcr_rd_deg$cluster]
  cellstate_pcr_rd_deg <- cellstate_pcr_rd_deg %>%
    dplyr::filter(!cluster %in% "Unresolved") %>%
    tidyr::unite(label, pCR_status, cluster) %>%
    dplyr::select(label, gene) %>%
    ruok::deframe_to_list()

  # cellstate_pcr_rd_deg <- ruok::deframe_to_list(cellstate_pcr_rd_deg[, c('cluster', 'gene')])
  str(cellstate_pcr_rd_deg)
  length(cellstate_pcr_rd_deg)
  class(cellstate_pcr_rd_deg)
  # table( names(module_content_use) %in% names(cellstate_pcr_rd_deg) )
  # idx <- ! names(cellstate_pcr_rd_deg) %in% c(names(module_content_use), 'Unresolved')
  # idx <- ! names(cellstate_pcr_rd_deg) %in% c('Unresolved')
  # cellstate_pcr_rd_deg <- cellstate_pcr_rd_deg[idx]
  names(module_content_use) <- paste0("MM_", names(module_content_use))
  names(cellstate_pcr_rd_deg)
  cellstate_pcr_rd_deg <- lapply(cellstate_pcr_rd_deg, head, 30)
  cellstate_pcr_rd_deg_MM <- cellstate_pcr_rd_deg[c(25:35, 77:87)]
  cellstate_pcr_rd_deg_TME <- cellstate_pcr_rd_deg[setdiff(names(cellstate_pcr_rd_deg), names(cellstate_pcr_rd_deg_MM))]
  stopifnot(length(unlist(cellstate_pcr_rd_deg)) == length(unlist(cellstate_pcr_rd_deg_MM)) + length(unlist(cellstate_pcr_rd_deg_TME)))

  cellstate_marker <- read_rds("~/project/tnbc_pre_atlas/rds_rna-integrate/pat102/atlas/findallmarker.ecotrait_feature/SigDEG_Pos_ecotrait_feature.df.rds")
  cellstate_marker$cluster <- dict_update_cellstatenames[as.character(cellstate_marker$cluster)]
  cellstate_marker <- cellstate_marker %>%
    dplyr::select(cluster, gene) %>%
    ruok::deframe_to_list()
  names(cellstate_marker) <- paste0("CS.", names(cellstate_marker))
  cellstate_marker <- lapply(cellstate_marker, head, 30)
  str(unlist(cellstate_marker))
  cellstate_marker_MM <- cellstate_marker[1:12]
  cellstate_marker_TME <- cellstate_marker[13:length(cellstate_marker)]
  stopifnot(length(unlist(cellstate_marker)) == length(unlist(cellstate_marker_MM)) + length(unlist(cellstate_marker_TME)))

  db <- c(module_content_use, nmf_genes, cellstate_pcr_rd_deg, cellstate_marker)
  str(unique(unlist(module_content_use)))
  str(unique(unlist(nmf_genes)))
  str(unique(unlist(cellstate_pcr_rd_deg)))
  str(unique(unlist(cellstate_marker)))
  str(unique(unlist(cellstate_marker_MM)))
  str(unique(unlist(cellstate_marker_TME)))
  str(unique(unlist(cellstate_pcr_rd_deg_MM)))
  str(unique(unlist(cellstate_pcr_rd_deg_TME)))

  str(db)
  length(db)
  str(unique(unlist(db)))
  names(db)
  # dir_res <- file.path(dir_proj, 'cancer_mm_archetype_genes_cellstate_pCRrelated')
  dir_res <- file.path(dir_proj, "cancer_mm_archetype_markers_CellstatePCRrelated_CellstateMarkers")
  fs::dir_create(dir_res)

  db_content <- list(
    MM = unique(unlist(module_content_use)),
    ARC = unique(unlist(nmf_genes)),
    CS = unique(unlist(cellstate_marker)),
    ResponseDEGs = unique(unlist(cellstate_pcr_rd_deg))
  )
  sapply(db_content, length)
  # MM          ARC           CS ResponseDEGs
  # 243          118          643          675
  library(grDevices)
  library(colorspace)
  p <- VennDiagram::venn.diagram(
    db_content,
    filename = NULL,
    col = rainbow_hcl(length(db_content)), cat.col = rainbow_hcl(length(db_content)),
    disable.logging = T, cat.default.pos = "text"
  )
  grid.newpage()
  pdf(file.path(dir_res, "vennplot.db_content.pdf"), width = 5, height = 5, useDingbats = F)
  grid.draw(p)
  dev.off()
}

p <- enframe(sapply(db, length)) %>%
  ggbarplot(y = "name", x = "value", fill = "black", color = NA) +
  labs(x = sprintf("nGene (%s)", length(unique(unlist(db))))) + theme(axis.text.y = element_text(size = 4))
p
ggsave(file.path(dir_res, "barplot.input_gene_num_per_features.pdf"), p, width = 4, height = 9)
#------ run univariate logistic regression ------
colnames(df_surv)
table(df_surv$pCR_status)

y_choice_str <- "pCR_status"
df_logistic <- df_surv[df_surv[, y_choice_str] %in% c("pCR", "RD"), ]
mat_logistic <- mat[, rownames(df_logistic)]

stopifnot(all.equal(df_logistic$patient, colnames(mat_logistic)))

tosearch <- unique(unlist(db))
cat(length(tosearch), " genes are asked to be searched...\n")
str(tosearch)
tmp <- setdiff(tosearch, rownames(mat_logistic))
cat(length(tmp), " genes are not present in training cohort so genes are removed...\n")
str(tmp)
mat_logistic <- mat_logistic[intersect(rownames(mat_logistic), tosearch), ]
cat(nrow(mat_logistic), " genes are used to search...\n")
str(rownames(mat_logistic))

rownames(mat_logistic) <- make.names(rownames(mat_logistic))
db_use <- rownames(mat_logistic)
table(duplicated(db_use))
df_logistic <- cbind(df_logistic, as.data.frame(t(mat_logistic)))

#------ run multivariable logistic regression ------
covariates <- c(db_use)
df_logistic$outcome_y <- df_logistic[, y_choice_str]
df_logistic$outcome_y <- forcats::fct_drop(df_logistic$outcome_y)
df_logistic$outcome_y <- as.numeric(df_logistic$outcome_y == "RD")

univ_formulas <- sapply(
  covariates,
  function(x) as.formula(paste("outcome_y ~ ", x))
)
univ_models <- lapply(univ_formulas, function(x) {
  glm(x, data = df_logistic, family = binomial)
})

univ_results <- lapply(
  univ_models,
  function(x) {
    x <- summary(x)
    return(x$coef[2, ])
  }
)
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
colnames(res) <- c("estimate", "std_error", "z", "pval")
res <- rownames_to_column(res, "univariate")
head(res)
res$univariate_source <- sapply(res$univariate, function(g) {
  out <- c()
  for (l in names(db)) {
    if (g %in% db[[l]]) {
      out <- c(out, l)
    }
  }
  return(paste0(sort(unique(out)), collapse = ","))
})

write_csv(res, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.report.csv", y_choice_str)
))


o <- res %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::arrange(estimate)
o <- res %>%
  dplyr::filter(pval < 0.01) %>%
  dplyr::arrange(estimate)
# o2 <- res %>% dplyr::filter(qval < 0.05) %>% dplyr::arrange(estimate)
# str(o2$univariate)

o$univariate <- str_replace(o$univariate, "\\.", "-") # HLA.DRA back to HLA-DRA

cat("gene number reduced from ", nrow(mat_logistic), "to", nrow(o), "\n")

# view(o)


write_csv(o, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.hits.csv", y_choice_str)
))


o_g <- o$univariate[o$estimate < 0]
o_g <- setdiff(o_g, "NPI")
str(o_g)
o_b <- o$univariate[o$estimate > 0]
o_b <- setdiff(o_b, "NPI")
str(o_b)
write_lines(o_g, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.hits_good_outcome.txt", y_choice_str)
))
write_lines(o_b, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.hits_bad_outcome.txt", y_choice_str)
))
o_combo <- list("GoodProg" = o_g, "BadProg" = o_b)
write_rds(o_combo, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.hits.rds", y_choice_str)
))
print(o_combo)


#------ Multivariable logistic regression ------
covariates <- make.names(unlist(o_combo))
str(covariates)
# if (F) {
if (length(covariates) > ncol(mat_logistic)) {
  cat("[warn] #feature > #sample!!")
  abs_odd_ratio_cutoff <- 2
  table(exp(abs(o$estimate)) >= abs_odd_ratio_cutoff)
  o_combo <- o %>%
    dplyr::filter(exp(abs(estimate)) >= abs_odd_ratio_cutoff) %>%
    dplyr::mutate(outcome_x = ifelse(estimate > 0, "BadProg", "GoodProg")) %>%
    dplyr::select(outcome_x, univariate) %>%
    ruok::deframe_to_list()
  print(o_combo)
  covariates <- make.names(unlist(o_combo))
  str(covariates)
}
all(covariates %in% colnames(df_logistic))
mulv_model <- NULL
mulv_model <- glm(
  as.formula(paste0("outcome_y ~ ", paste0(covariates, collapse = "+"))),
  data = df_logistic, family = binomial
)
mulv_model_res <- as.data.frame(summary(mulv_model)$coef)
mulv_model_res <- rownames_to_column(mulv_model_res, "variate")
mulv_model_res <- mulv_model_res %>% dplyr::arrange(desc(Estimate))

write_rds(mulv_model_res, file.path(
  dir_res,
  sprintf("%s.multivariable_logistic_regression_model_summary.rds", y_choice_str)
))
write_csv(mulv_model_res, file.path(
  dir_res,
  sprintf("%s.multivariable_logistic_regression_model_summary.csv", y_choice_str)
))
mulv_model_genes <- sort(covariates)
write_lines(mulv_model_genes, file.path(
  dir_res,
  sprintf("%s.multivariable_logistic_regression_model_genes.txt", y_choice_str)
))
# mulv_model_res[mulv_model_res$`Pr(>|z|)` < 0.05, ]

write_rds(mulv_model, file.path(
  dir_res,
  sprintf("%s.multivariable_logistic_regression_model.rds", y_choice_str)
))

# view(mulv_model_res)
gene_source <- structure(o$univariate_source, names = o$univariate)
mulv_model_res$source <- gene_source[mulv_model_res$variate]
mulv_model_res_simple <- mulv_model_res %>% dplyr::select(variate, Estimate, source)

write_csv(mulv_model_res_simple, file.path(
  dir_res,
  sprintf("%s.multivariable_logistic_regression_model_deliver.csv", y_choice_str)
))
o$used_in_model <- o$univariate %in% covariates
res$used_in_model <- res$univariate %in% covariates
write_csv(o, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.hits.csv", y_choice_str)
))
write_csv(res, file.path(
  dir_res,
  sprintf("%s.Univariate_logistic_regression.report.csv", y_choice_str)
))
colnames(mulv_model_res) <- c("variate", "Estimate", "SDE", "z", "Pr", "source")
p <- mulv_model_res %>%
  dplyr::filter(!str_detect(variate, "Intercept")) %>%
  dplyr::arrange(desc(Estimate)) %>%
  ggplot(., aes(y = reorder(variate, Estimate), x = Estimate)) +
  geom_errorbar(aes(ymin = Estimate - 2 * SDE, ymax = Estimate + 2 * SDE)) +
  geom_segment(aes(
    y = reorder(variate, Estimate), x = 0,
    yend = reorder(variate, Estimate), xend = Estimate
  )) +
  geom_point(size = 3, aes(color = ifelse(Estimate > 0, "RD", "pCR"))) +
  scale_color_manual(values = c(
    "pCR" = "#53D43F",
    "RD" = "#811C9A"
  )) +
  geom_vline(xintercept = 0) +
  labs(x = "beta coefficient", y = "gene")
p <- p + theme_pubr() + rremove("legend")
p
ggsave(
  file.path(
    dir_res,
    sprintf("%s.multivariable_logistic_regression_model_deliver.loliplot.pdf", y_choice_str)
  ),
  p,
  width = 4, height = 5, useDingbats = F
)

#------ performance on its own training data ------
pred_prob <- predict(mulv_model, df_logistic[, covariates], type = "response")
head(pred_prob)
hist(pred_prob)

library(ggpubr)
p <- data.frame(prob = pred_prob, obs = df_logistic$outcome_y) %>%
  ggboxplot(data = ., x = "obs", y = "prob") +
  stat_compare_means()
pdf(file.path(
  dir_res,
  sprintf("%s.multivariable_logistic_regression_model.pred_prob_train.pdf", y_choice_str)
), useDingbats = F)
print(p)
dev.off()

