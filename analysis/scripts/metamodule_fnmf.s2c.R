# Propose metaprograms by various choices of hclust k and linkage
#  
# Yun Yan (yun.yan@uth.tmc.edu)
#
suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(fs)
    library(fastcluster)
    library(dendextend); library(dendsort)
    source('util.nmf.viz.R')
})

dir_proj <- file.path(
    '.', 
    'metamodule_fnmf', 
    'all_ranks')

z_opts <- c('jaccard')
z <- z_opts[1]
dir_res <- file.path(dir_proj, z, 'hclust'); fs::dir_create(dir_res)
mat <- read_rds( file.path(dir_proj, sprintf('%s_mat.all_ranks.in_use.rds', z) ))

print(dim(mat))

if (FALSE) {
    ## TESTING CODE MODE
    ## randomly subset 1000 factors
    idx <- sample(1:nrow(mat), size = 1000, replace = F)
    mat <- mat[idx, idx]
}

#------ nmf ~ metamodules by hclust ------
## Try all linkage. Ref: https://stats.stackexchange.com/a/217742
hclust_obj <- NULL
hclust_linkage_opt <- c('complete', 'ward.D', 'average', 'ward.D2')
n_modules_opt      <- c(12, 14, 16, 18, 20, 24, 26, 30)

for (hclust_method in hclust_linkage_opt) {
    cat(hclust_method, '...')
    f_nmf_hclust <- file.path(dir_res, sprintf('NMFs_%s_hclust.%s.rds', z, hclust_method))
    if (!file.exists(f_nmf_hclust)) {
        hclust_obj <- try( fastcluster::hclust(
            d=as.dist(1-mat), method=hclust_method) ) 
        if ('try-error' %in% class(hclust_obj)) {next()}
        write_rds(hclust_obj, f_nmf_hclust)
    } else {
        hclust_obj <- read_rds(f_nmf_hclust)
    }
    for (n_modules in n_modules_opt) {
        cat('considering', n_modules, 'meta-modules...')
        # plot(color_branches(hclust_obj, k=n_modules, groupLabels = T))
        hclust_cnmf_mem <- dendextend::cutree(
            hclust_obj, k = n_modules, 
            order_clusters_as_data=FALSE)
        hclust_cnmf_mem <- structure(
            paste0('M', str_pad(hclust_cnmf_mem, width = 2, pad = '0')),
            names = names(hclust_cnmf_mem))
        module_names <- sort(unique(hclust_cnmf_mem))
        
        
        p2 <- heatmap_cor_mat(
            cor_mat = mat[names(hclust_cnmf_mem), names(hclust_cnmf_mem)],
            hclust_obj = NULL,
            # col = heatmap_color_correlation_fun,
            col = colorRamp2(
                seq(from=0, to=1, length.out = 5),
                colors = sequential_hcl(n=5, palette = 'Oslo', rev = F)
            ),
            cor_name = z, run_name = 'Modules defined by hclust',
            column_gap = unit(0, 'mm'), row_gap = unit(0, 'mm'), border=T,
            row_title_rot = 0, column_title_rot=90,
            row_split=hclust_cnmf_mem,
            column_split=hclust_cnmf_mem,
            raster_by_magick = T,
            topdf = file.path(dir_res,
                              sprintf('heatmap.NMFs_%s_to_modules.M%s.by_hclust_%s.pdf',
                                      z, n_modules, hclust_method)),
            use_raster = T, raster_device = 'CairoPNG')
        stopifnot(all(names(hclust_cnmf_mem) %in% colnames(mat)))
        
        ## How many factors per mema-module?
        pdf(file.path(dir_res, sprintf('barplot.%s.num_factors_per_module.M%s.hclust_%s.pdf', 
                                       z, n_modules, hclust_method)), 
            width = 7, height = 4)
        barplot(table(hclust_cnmf_mem))
        dev.off()
        ## How many samples per meta-module?
        hclust_cnmf_mem <- enframe(hclust_cnmf_mem, name = 'index_name', value = 'meta_module')
        hclust_cnmf_mem <- hclust_cnmf_mem %>% tidyr::separate(index_name, c('sample', 'Rank', 'NMF'), sep='_', remove=F)
        head(hclust_cnmf_mem)
        df <- hclust_cnmf_mem %>% dplyr::group_by(meta_module) %>%
            dplyr::summarise(n=length(unique(sample))) %>% dplyr::ungroup()
        p <- ggplot(df, aes(x=meta_module, y=n)) + geom_bar(stat = 'identity') +
            geom_text(aes(label=n)) + ggpubr::rotate_x_text(90) +
            labs(y='number of sample', caption = sprintf('avg %.0f samples per meta module', mean(df$n)))
        ggsave(file.path(dir_res, sprintf('barplot.%s.num_samples_per_module.M%s.hclust_%s.pdf', 
                                          z, n_modules, hclust_method)), 
               p, width = 7, height = 4)
        ## Export membership per meta-module
        write_csv(hclust_cnmf_mem, 
                  file.path(dir_res, 
                            sprintf('metamodule_membership.%s.%s.hclust_%s.csv', 
                                    z, n_modules, hclust_method)))
        n_modules
    }
    hclust_obj <- NULL
    cat('[done]\n')
}; cat('DONE\n')
# parallel::stopCluster()
#------ survey the choice of number of meta modules ------
# ...
    