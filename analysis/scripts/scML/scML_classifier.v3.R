#--------------------------
# Build a classifier using the meta-module usage and TME cell state frequencies
# purely based on logistic regression
#--------------------------
library(tidyverse)
library(pROC)
library(caret)
library(LogicReg)
library(ComplexHeatmap); library(magrittr)
set.seed(2022)

#--------------------------
# General training settings 
#--------------------------
# 10-fold CV resampling with subsampling via ROSE approach
trCtrol <- trainControl(method = 'cv', number = 10, 
                        sampling = 'rose',
                        classProbs = TRUE, summaryFunction = twoClassSummary)
# 3-fold CV resampling without subsampling
trCtrol2 <- trainControl(method = 'cv', number = 3,
                         classProbs = TRUE, summaryFunction = twoClassSummary)
# 5-times repeated 5-fold cv 
trCtrolRepCV <- trainControl(method = 'repeatedcv', 
                             number = 5, repeats = 5, 
                             sampling = 'rose', 
                             classProbs = TRUE, summaryFunction = twoClassSummary)

#--------------------------
# Prepare data compatible with caret package
#--------------------------
if (T) {

  df_long <- read_rds(
    file.path(
      '~/project/tnbc_pre_atlas/rds_rna-integrate/',
      'pat102/lv01.aneuploidy_tri_type.aneuploid.pure5',
      'viz_signature_MM_alt_clean_byscore_wardD2',
      'hybrid',
      'df_patient_cellpct.rds'))
  mat <- read_rds(
    file.path(
      '~/project/tnbc_pre_atlas/rds_rna-integrate/',
      'pat102/lv01.aneuploidy_tri_type.aneuploid.pure5',
      'viz_signature_MM_alt_clean_byscore_wardD2',
      'hybrid',
      'mat_patient_cellpct.rds'))
  
}

if (T) {
  setwd('/volumes/USR1/yyan/project/tnbc_pre_atlas')
  
  dir_proj <-file.path(
    '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102',
    'ecotype', 'use_tumor_hybrid_frac', 'allpatient_with_archetypes')
  dir_res <- file.path(dir_proj, 'scML_classifier_v3')
  fs::dir_create(dir_res)
  
  df_long$cellstate <- make.names(df_long$cellstate)
  colnames(mat) <- make.names(colnames(mat))
}


df_wide <- as.data.frame(mat)
head(df_wide)
pcr_patient_info <- deframe( unique( df_long[, c('patient', 'pCR_status')] ) )
nmf_patient_info <- deframe( unique( df_long[, c('patient', 'archetype')] ) )
print(table(pcr_patient_info))
head(pcr_patient_info)
head(nmf_patient_info)
table(nmf_patient_info, useNA='ifany')





vec_to_indicator_df <- function(x){
  x_names <- names(x)
  y_lvs <- gtools::mixedsort( unique(x) )
  out <- sapply(y_lvs, function(lv) {x == lv})
  out <- as.data.frame(out)
  colnames(out) <- y_lvs
  rownames(out) <- x_names
  return(out)
  
}
nmf_indicator <- vec_to_indicator_df(nmf_patient_info)
nmf_indicator <- nmf_indicator * 1

df <- df_wide

#---------------------------
# Prepare inputs              ----    
#---------------------------
if (T) {
  # Take care of NA values
  # remove features having a lot NAs
  # fill NA with average
  n_pat <- nrow(mat); n_var <- ncol(mat)
  n_na_per_feature <- apply(mat, 2, function(v) sum(is.na(v)))
  feature_nonna <- colnames(mat)[n_na_per_feature <= 10]
  feature_nonna
  
  df <- mat[, feature_nonna]
  df <- apply(df, 2, function(v) {
    v_avg <- mean(v, na.rm=T)
    v[is.na(v)] <- v_avg
    v
  })
  df <- as.data.frame(df)
}
if (F) {
  # if add the indicator variable
  df <- cbind(df, nmf_indicator[rownames(df), ])
}
df$Y <- pcr_patient_info[rownames(df_wide)]
print(table(df$Y))
df <- df[df$Y %in% c('pCR', 'RD'), ]
df$Y <- factor(as.character(df$Y), levels = c('pCR', 'RD'))

set.seed(2022)
pct_for_training <- 0.7
idx <- sample(seq_len(nrow(df)), 
              size=round(nrow(df) * pct_for_training), replace = F)
idx <- seq_len(nrow(df)) %in% idx
table(idx)
dim(df)
learn <- df[idx, ]
valid <- df[!idx, ]
head(learn) ## dataset to learn (internal training and testing)
head(valid) ## dataset to validate (pure AUC evaluation)

covariates <- setdiff( colnames(df), c('Y', 'ARC1', 'ARC2', 'ARC3', 'ARC4' ))
univ_formulas <- sapply(
  covariates,
  function(x) as.formula(paste('Y ~ ', x)))
univ_models <- lapply( univ_formulas, function(x){
  glm(x, data = df, family = binomial)}) 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         return( x$coef[2, ] )
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
colnames(res) <- c('estimate', 'std_error', 'z', 'pval')
res <- rownames_to_column(res, 'univariate')
head(res); nrow(res)
res$is_sig <- res$pval < 0.05
view(res)
features <- res$univariate[res$pval < 0.1]
res$direction <- ifelse(res$estimate > 0, 'RD', 'pCR')
str(features)

res$used_for_building_model <- ifelse(res$univariate %in% features, 'yes', 'no')
res <- res %>% dplyr::arrange(pval)
write_csv(res, file.path(dir_res, 'logistic_univariate_beta_coefficient_on_all_dataset.csv'))


theme_set(theme_pubr(legend = 'right'))
p <- ggplot(res, aes(y=reorder(univariate, estimate), x=estimate)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(y = reorder(univariate, estimate), x = 0, 
                   yend = reorder(univariate, estimate), xend = estimate)) +
  geom_point(aes(color=direction, shape=is_sig), 
             size = 2.5) +
  scale_color_manual(values=c('pCR'='#53D43F','RD'='#811C9A')) +
  # scale_color_manual(values=c(`TRUE`='magenta',`FALSE`='white')) +
  labs(x = 'beta coefficient', y='cell state')
ggsave(file.path(dir_res, 'logistic_univariate_beta_coefficient_on_all_dataset.pdf'), 
       p, width = 5, height = 7, useDingbats = F)

#--------------------------
# Training ----
#--------------------------

print(table(learn$Y))
do_train <- T
seed_num <- 22
seed_num <- 42
if (do_train) {
  set.seed(seed_num)
  metric <- "ROC"
  mdl_log <- train(as.formula(paste0('Y~', paste(features, collapse = '+'))), 
                   data = learn, 
                   method = "glm",
                   family = 'binomial', 
                   trControl = trCtrolRepCV, 
                   # preProcess = c('center', 'scale'),
                   metric = metric)
  set.seed(seed_num)
  mdl_lda <- train(as.formula(paste0('Y~', paste(features, collapse = '+'))),
                   data = learn, 
                   method = "lda",  
                   trControl = trCtrolRepCV, 
                   # preProcess = c('center', 'scale'),
                   metric = metric)
  set.seed(seed_num)
  rfGrid <- expand.grid(mty=c(3, 10, 20))
  mdl_rf <- train(as.formula(paste0('Y~', paste(features, collapse = '+'))),
                  data = learn, 
                  method = "rf",  
                  trControl = trCtrolRepCV, 
                  # preProcess = c('center', 'scale'),
                  # tuneGrid = rfGrid,
                  metric = metric,
                  verbose = TRUE)
  
  write_rds(mdl_log, file = file.path(dir_res, 'mdl_log.rds'))
  write_rds(mdl_lda, file = file.path(dir_res, 'mdl_lda.rds'))
  write_rds(mdl_rf, file = file.path(dir_res, 'mdl_rf.rds'))
} else {
  mdl_log <- read_rds(file = file.path(dir_res, 'mdl_log.rds'))
  mdl_lda <- read_rds(file = file.path(dir_res, 'mdl_lda.rds'))
  mdl_rf  <- read_rds(file = file.path(dir_res, 'mdl_rf.rds'))
}    
#--------------------------
# Validation and visualization ----
#--------------------------
report_valid_roc <- function(model, dt, dt_lab_fct) {
  class_focus <- levels(dt_lab_fct)[2]
  roc_obj <- roc(dt_lab_fct, 
                 predict(model, dt, type = 'prob')[, class_focus])
  ci(roc_obj)
}
valid_roc_obj <- function(model, dt, dt_lab_fct) {
  class_focus <- levels(dt_lab_fct)[2]
  roc_obj <- roc(dt_lab_fct, 
                 predict(model, dt, type = 'prob')[, class_focus])
  roc_obj
}
plot_models_valid_roc <- function(model_list, dt, dt_lab_fct, main = NULL) {
  n <- length(model_list)
  mycolors <- RColorBrewer::brewer.pal(ifelse(n < 3, 3, n), 'Set1')
  i <- 1
  model_auc <- rep(0, n)
  for (model in model_list) {
    message(names(model_list)[[i]])
    
    roc_obj <- valid_roc_obj(model, dt, dt_lab_fct)
    if (i == 1) {
      plot.roc(roc_obj, col = mycolors[i], main = main, asp = 1)
    } else {
      lines.roc(roc_obj, col = mycolors[i])
    }
    model_auc[i] <- as.numeric(roc_obj$auc)
    i <- i + 1
  }
  model_auc <- format(model_auc, digits = 3)
  legend_labels <- paste0(names(model_list), " (", model_auc, ")")
  legend("bottomright", legend = legend_labels, col = mycolors[1:n], lwd = 2)
  invisible(0)
}

mdl_zoo <- list('Logreg' = mdl_log,
                'LDA'    = mdl_lda,
                'RF'     = mdl_rf)

pdf(file.path(dir_res, 'roc_zoo.validation.pdf'), width = 6, height = 6, useDingbats = F)
plot_models_valid_roc(mdl_zoo, dt=valid[, features], dt_lab_fct = valid$Y)
dev.off()

pdf(file.path(dir_res, 'roc_zoo.training.pdf'), width = 6, height = 6, useDingbats = F)
plot_models_valid_roc(mdl_zoo, dt=learn[, features], dt_lab_fct = learn$Y)
dev.off()

library(ComplexHeatmap)
for (i in seq_along(mdl_zoo)){
  mdl <- mdl_zoo[[i]]
  pre_prob <- predict(mdl, newdata = valid[, features], type = "prob")
  pre_cat <- predict(mdl, newdata = valid[, features])
  table(pre_prob$pCR > 0.5)
  reprt <- table(unname(pre_cat), valid$Y)
  reprt <- as.data.frame.matrix(reprt)
  
  pdf(file.path(
    dir_res, 
    sprintf('boxplot_valid_dataset.%s.pdf', names(mdl_zoo)[[i]])),
    width = 4, height = 4, onefile = F, useDingbats = F)  
  boxplot(pre_prob$RD~valid$Y, main=names(mdl_zoo)[[i]],ylab='risk score') ## this is better
  dev.off()
  
  pdf(file.path(
    dir_res, 
    sprintf('contincy_table_valid_dataset.%s.pdf', names(mdl_zoo)[[i]])),
    width = 2.5, height = 2.5, onefile = F, useDingbats = F)
  draw(
    Heatmap(reprt, cluster_rows = F, cluster_columns = F, 
            heatmap_width = unit(2, 'inch'),
            heatmap_height = unit(2, 'inch'),
            cell_fun = function(j, i, x, y, width, height, fill) {
              
              grid.text(reprt[i, j], x, y, 
                        gp = gpar(col='black'))
            }), show_heatmap_legend=F)
  dev.off()
  
}

for (i in seq_along(mdl_zoo)){
  mdl <- mdl_zoo[[i]]
  pre_prob <- predict(mdl, newdata = learn[, features], type = "prob")
  pre_cat <- predict(mdl, newdata = learn[, features])
  table(pre_prob$pCR > 0.5)
  reprt <- table(unname(pre_cat), learn$Y)
  reprt <- as.data.frame.matrix(reprt)
  
  pdf(file.path(
    dir_res, 
    sprintf('boxplot_learn_dataset.%s.pdf', names(mdl_zoo)[[i]])),
    width = 4, height = 4, onefile = F, useDingbats = F)  
  boxplot(pre_prob$RD~learn$Y, main=names(mdl_zoo)[[i]],ylab='risk score') ## this is better
  dev.off()
  
  pdf(file.path(
    dir_res, 
    sprintf('contincy_table_learn_dataset.%s.pdf', names(mdl_zoo)[[i]])),
    width = 2.5, height = 2.5, onefile = F, useDingbats = F)
  draw(
    Heatmap(reprt, cluster_rows = F, cluster_columns = F, 
            heatmap_width = unit(2, 'inch'),
            heatmap_height = unit(2, 'inch'),
            cell_fun = function(j, i, x, y, width, height, fill) {
              
              grid.text(reprt[i, j], x, y, 
                        gp = gpar(col='black'))
            }), show_heatmap_legend=F)
  dev.off()
  
}
pal_pcr <- c('pCR'='#53D43F','RD'='#811C9A')
for (i in seq_along(mdl_zoo)){
  mdl <- mdl_zoo[[i]]
  # var_importance <- varImp(mdl, scale = T)
  var_importance <- varImp(mdl, scale = F)

  var_direction <- t( sapply(features, function(x) {
    tapply(df[, x], df$Y, mean)}) )
  var_direction <- var_direction[, 'RD'] - var_direction[, 'pCR']
  var_direction <- factor(ifelse(var_direction > 0, 'RD', 'pCR'), levels=c('pCR', 'RD'))
  
  pdf(file.path(
    dir_res, 
    sprintf('feature_imporantce.%s.pdf', names(mdl_zoo)[[i]])), 
    width = 5, height = 6, onefile = F, useDingbats = F)
  
  # print(plot(var_importance, top = pmin(length(features), 20), cex=1))
  
  p <- ggplot(var_importance, top = pmin(length(features), 20))
  var_importance <- p$data
  var_importance$direction <- var_direction[as.character(var_importance$Feature)]
  
  p <- var_importance %>% dplyr::top_n(., n=pmin(length(features), 20), wt=Importance) %>%
    ggplot(aes(y=Feature, x=Importance)) +
    geom_vline(xintercept = 0) +
    geom_segment(aes(y = Feature, x = 0, 
                     yend = Feature, xend = Importance)) +
    geom_point(aes(color=direction), size=3) +
    labs(x='feature importance', y='cell state or meta-trait') +
    scale_color_manual(values=pal_pcr) + rremove('legend')
  print(p)
  dev.off()
}
