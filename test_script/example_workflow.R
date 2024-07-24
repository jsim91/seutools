# rm(list = ls()); gc()

# example workflow
library(seutools)
library(Seurat)
library(Matrix)

seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")

und_cl <- paste0("Undecided",1:5)

# testing seutools::seurat_mast()
if(F) {
  und_seu_mast <- seurat_mast(seurat_object = seu, freq_expressed = 0.1, fc_threshold = log2(1.5),
                              test_per_cluster = TRUE, test_clusters = und_cl, cluster_column = "cell_type",
                              test_per_category = FALSE, test_categories = c("younger","older"),
                              category_column = "age_group", test_per_condition = FALSE,
                              test_condition = "all", condition_column = "condition", pid_column = "pid")
}

# c('ISG_Naive_CD4','IFN_CM_CD4','ISG_EM_CD4','ISG NK','ISG_Naive_CD8','ISG_CTL_CD8','ISG_Mono')
if(F) {
  seu_wilc_select <- seurat_dge(seurat_object = seu, dge_method = "wilcox", freq_expressed = 0.1,
                                fc_threshold = log2(1.5), test_per_cluster = TRUE,
                                test_clusters = c('ISG_Naive_CD4','IFN_CM_CD4','ISG_EM_CD4','ISG NK',
                                                  'ISG_Naive_CD8','ISG_CTL_CD8','ISG_Mono'),
                                cluster_column = "cell_type", test_per_category = FALSE,
                                test_per_condition = FALSE, test_condition = "all")
  # saveRDS(object = seu_wilc_select, file = "J:/U54_grant/sc/inputs/select_wilcox_dge.rds")
} else {
  seu_wilc_select <- readRDS("J:/U54_grant/sc/inputs/select_wilcox_dge.rds")
}
seu_wilc_collapsed <- do.call(rbind, seu_wilc_select[[1]])
vol1 <- seutools:::plot_volcano(dge_input = seu_wilc_collapsed, plot_clusters = "all",
                                gene_set = NA, prio_top_genes = 0, pval_threshold = 1,
                                table_height = 50, fc_threshold = log2(1.5),
                                de_method = "seurat_presto")


# testing seutools::seurat_feature_overlay()
library(Matrix)
library(seutools)
library(Seurat)
library(SeuratObject)

seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")

sof_adt <- seutools::seurat_feature_overlay(seurat_object = seu_adt, text_expansion = 1.5, draw_legend = TRUE)
sof_adt_tiled <- seutools:::tile_plots(plotlist = sof_adt, n_row = 2, n_col = 2, rm_legend = FALSE)

length(sof_adt_tiled)

pdf(file = "J:/U54_grant/sc/out_figures/flu_adt_umap_dsb_normed.pdf", width = 16, height = 16)
lapply(X = sof_adt_tiled, FUN = function(x) x)
dev.off()

sof_gdt <- seutools::seurat_feature_overlay(seurat_object = seu_adt, text_expansion = 1.5,
                                            draw_legend = FALSE, plot_features = c("Hu.TCR.Vd2","Hu.CD1c"),
                                            label_clusters = "all", text_expansion_annotate = 0.5)
pdf(file = "J:/U54_grant/sc/out_figures/flu_adt_labeled_umap.pdf", width = 14, height = 14)
sof_gdt$Hu.TCR.Vd2
dev.off()


# testing seutools::seurat_test_clusters()
library(Seurat)
library(Matrix)

seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")

stc_compares <- list(c('media\nyounger','stim\nyounger'),
                     c('media\nolder','stim\nolder'))

stc4 <- seurat_test_clusters(seurat_object = seu, test_by_column = "condition_age", pid_column = "pid",
                             cluster_column = "cell_type", cell_barcode_column = "barcode",
                             order_x = c('media\nyounger','stim\nyounger','media\nolder','stim\nolder'),
                             bin_colors = c("azure4", "seagreen", "azure4", "royalblue"),
                             subtract_background = FALSE, comparison_list = stc_compares,
                             backround_condition = "media", y_axis_subset = "PBMC", connect_points = TRUE,
                             data_paired = TRUE, return_plots = TRUE, return_plot_data = FALSE,
                             coord_stretch_factor = 0.11, text_size_factor = 1, shape_key = NULL)
stc4_tiled <- seutools:::tile_plots(plotlist = stc4, n_row = 2, n_col = 2, rm_legend = FALSE)
pdf(file = "J:/U54_grant/sc/out_figures/flu_test_condition_wide.pdf", width = 16, height = 16)
lapply(X = stc4_tiled, FUN = function(x) x)
dev.off()

my_colors <- c('IFN_CM_CD4' = '#96CA2D',
               'ISG NK' = '#E53935',
               'ISG_CTL_CD8' = '#FF6517',
               'ISG_EM_CD4' = '#9250BC',
               'ISG_Naive_CD4' = '#F2B705',
               'ISG_Naive_CD8' = '#04BFBF',
               'ISG_Mono' = '#607D8B')

stc5 <- seurat_test_clusters(seurat_object = seu, test_by_column = "condition_age", pid_column = "pid",
                             cluster_column = "cell_type", cell_barcode_column = "barcode",
                             order_x = c('media\nyounger','stim\nyounger','media\nolder','stim\nolder'),
                             bin_colors = c("azure4", "seagreen", "azure4", "royalblue"),
                             subtract_background = FALSE, comparison_list = stc_compares, y_as_log = TRUE,
                             backround_condition = "media", y_axis_subset = "PBMC", connect_points = TRUE,
                             data_paired = TRUE, return_plots = TRUE, return_plot_data = FALSE,
                             coord_stretch_factor = 0.11, text_size_factor = 0.6, shape_key = NULL,
                             coordinate_title_color = my_colors)
reac_cl <- stc5[c('ISG_Naive_CD4','IFN_CM_CD4','ISG_EM_CD4','ISG NK','ISG_Naive_CD8','ISG_CTL_CD8','ISG_Mono')]
reac_cl_arr <- ggpubr::ggarrange(plotlist = reac_cl, nrow = 2, ncol = 4)
# reac_cl_arr

# testing seutools::seurat_tile_reduction()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")

tr <- seutools::seurat_tile_reduction(seurat_object = seu, condition_column = "condition", cluster_column = "cell_type", reduction = "umap",
                                      color_clusters = "all", label_clusters = "all",
                                      pt_alpha = 0.1, text_expansion = 1, pt_size = 1, color_seed = 123,
                                      postfix_title_string = NA,
                                      force_xlim = FALSE, force_ylim = FALSE, return_as_list = FALSE,
                                      plot_order = c(1,2), annotation_method = "repel",
                                      override_color_aes = NA, frameon = FALSE)
ggsave(filename = "flu_umap_by_condition.png", plot = tr, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 16, height = 8, units = "in", dpi = 600, limitsize = F, bg = "white")

seu_younger <- subset(x = seu, subset = age_group == "younger"); seu_younger <- AddMetaData(seu_younger, paste0(seu_younger$condition," (younger)"), "stim_age_group")
seu_older <- subset(x = seu, subset = age_group == "older"); seu_older <- AddMetaData(seu_older, paste0(seu_older$condition," (older)"), "stim_age_group")

tr_act <- seutools::seurat_tile_reduction(seurat_object = seu, condition_column = "condition", cluster_column = "cell_type", reduction = "umap",
                                          color_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                          label_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                          pt_alpha = 0.1, text_expansion = 1.5, pt_size = 1, color_seed = 123,
                                          postfix_title_string = NA, force_colors = my_colors,
                                          force_xlim = FALSE, force_ylim = FALSE, return_as_list = FALSE,
                                          plot_order = c(1,2), annotation_method = "repel",
                                          override_color_aes = NA, frameon = FALSE)
ggsave(filename = "flu_umap_by_condition_select_clusters.png", plot = tr_act, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 16, height = 8, units = "in", dpi = 600, limitsize = F, bg = "white")

tr_act_y <- seutools::seurat_tile_reduction(seurat_object = seu_younger, condition_column = "stim_age_group", cluster_column = "cell_type", reduction = "umap",
                                            color_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                            label_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                            pt_alpha = 0.1, text_expansion = 1.5, pt_size = 1, color_seed = 123,
                                            postfix_title_string = NA, force_colors = my_colors,
                                            force_xlim = FALSE, force_ylim = FALSE, return_as_list = FALSE,
                                            plot_order = c(1,2), annotation_method = "repel",
                                            override_color_aes = NA, frameon = FALSE)
tr_act_o <- seutools::seurat_tile_reduction(seurat_object = seu_older, condition_column = "stim_age_group", cluster_column = "cell_type", reduction = "umap",
                                            color_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                            label_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                            pt_alpha = 0.1, text_expansion = 1.5, pt_size = 1, color_seed = 123,
                                            postfix_title_string = NA, force_colors = my_colors,
                                            force_xlim = FALSE, force_ylim = FALSE, return_as_list = FALSE,
                                            plot_order = c(1,2), annotation_method = "repel",
                                            override_color_aes = NA, frameon = FALSE)
arr_act_yo <- ggpubr::ggarrange(plotlist = list(tr_act_y, tr_act_o), nrow = 2, ncol = 1)
ggsave(filename = "flu_umap_by_condition_select_clusters_younger_older.png", plot = arr_act_yo, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 16, height = 16, units = "in", dpi = 600, limitsize = F, bg = "white")

tr_act_pair <- seutools::seurat_tile_reduction(seurat_object = seu, condition_column = "condition", cluster_column = "cell_type", reduction = "umap",
                                               color_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                               label_clusters = c('IFN_CM_CD4','ISG_EM_CD4','ISG_Naive_CD4','ISG_CTL_CD8','ISG_Naive_CD8','ISG_Mono','ISG NK'),
                                               pt_alpha = 0.1, text_expansion = 1.5, pt_size = 1, color_seed = 123,
                                               postfix_title_string = NA, force_colors = my_colors,
                                               force_xlim = FALSE, force_ylim = FALSE, return_as_list = FALSE,
                                               plot_order = c(1,2), annotation_method = "repel",
                                               override_color_aes = NA, frameon = FALSE)
arr_fig <- ggpubr::ggarrange(plotlist = list(tr_act_pair, reac_cl_arr), nrow = 2, ncol = 1, heights = c(0.5,0.5))
ggsave(filename = "flu_umap_test_pair.png", plot = arr_fig, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 16, height = 16, units = "in", dpi = 600, limitsize = F, bg = "white")


# testing seutools::seurat_feature_violin()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")

feature_genes <- c("LAG3","IL32","IFIT1","IFITM1","CD8A","CD79A","NCAM1","CLEC12A")
sfv_genes <- seutools::seurat_feature_violin(seurat_object = seu, plot_features = feature_genes,
                                             categorical_column = "cell_type",
                                             plot_categorical_types = "all", assay = "RNA", text_expansion = 1,
                                             nudge_nonzero = 0.35, y_limit_expansion_factor = 0.5,
                                             condition = "media", condition_cat = "condition")
pdf(file = "J:/U54_grant/sc/out_figures/flu_violin_genes_1.pdf", width = 10, height = 12)
lapply(X = sfv_genes, FUN = function(x) x)
dev.off()

feature_prot <- c("Hu.CD4-RPA.T4","Hu.CD8","Hu.CD56","Hu.CD57","Hu.CD11c","Hu.CD62L","Hu.CLEC12A","Hu.CD123")
sfv_prot <- seutools::seurat_feature_violin(seurat_object = seu_adt, plot_features = feature_prot,
                                            categorical_column = "cell_type",
                                            plot_categorical_types = "all", assay = "ADT", text_expansion = 1,
                                            nudge_nonzero = 0.35, y_limit_expansion_factor = 0.5,
                                            condition = "media", condition_cat = "stim")
pdf(file = "J:/U54_grant/sc/out_figures/flu_violin_adt_1.pdf", width = 10, height = 12)
lapply(X = sfv_prot, FUN = function(x) x)
dev.off()


# testing seutools::seurat_footprint()
leiden_foot <- seutools::seurat_footprint(seurat_object = seu, cluster_column = "leiden", pid_column = "pid", condition_column = "condition",
                                          media_condition = "media", subtract_media = TRUE, color_by_column = "age_group",
                                          scale.factor = 1000, pca_fraction_variance = 0.95, umap_n_neighbors = 5, cluster_data = FALSE,
                                          leiden_resolution = 0.5, umap_min_dist = 0.25, report_values_as = "normalized counts",
                                          feature_reduction_method = "pca", reduction = "umap")
ggsave(filename = "flu_umap_leiden_cluster_footprint.png", plot = leiden_foot, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 24, height = 22, units = "in", dpi = 600, limitsize = F, bg = "white")

annotation_foot <- seutools::seurat_footprint(seurat_object = seu, cluster_column = "cell_type", pid_column = "pid", condition_column = "condition",
                                              media_condition = "media", subtract_media = TRUE, color_by_column = "age_group",
                                              scale.factor = 1000, pca_fraction_variance = 0.95, umap_n_neighbors = 5, cluster_data = FALSE,
                                              leiden_resolution = 0.5, umap_min_dist = 0.25, report_values_as = "normalized counts",
                                              feature_reduction_method = "pca", reduction = "umap")
ggsave(filename = "flu_umap_cell_type_footprint.png", plot = annotation_foot, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 24, height = 22, units = "in", dpi = 600, limitsize = F, bg = "white")

res1p5_by_group <- read.csv(file = "J:/10x/JOflu/scbp/all/all_by_type_res1p5.csv", check.names = FALSE, row.names = 1)
clnum <- res1p5_by_group$leiden_res1p5; names(clnum) <- res1p5_by_group$barcode
mapped_num <- clnum[seu@meta.data$barcode]
mapped_num[which(seu$annotated_type=="Undecided2")] <- "U2"
mapped_num[which(seu$annotated_type=="Undecided3")] <- "U3"
type_df <- data.frame(barcode = names(mapped_num), cluster = as.character(mapped_num))
write.csv(x = type_df, file = "J:/10x/JOflu/scbp/all/lineages_reclustered_res1p5.csv", row.names = FALSE)
seu <- SeuratObject::AddMetaData(object = seu, metadata = as.character(mapped_num), col.name = "res1p5_by_type")

annotation_type_1p5 <- seutools::seurat_footprint(seurat_object = seu, cluster_column = "res1p5_by_type", pid_column = "pid",
                                                  condition_column = "condition", media_condition = "media", subtract_media = TRUE,
                                                  color_by_column = "age_group", scale.factor = 1000, pca_fraction_variance = 0.95,
                                                  umap_n_neighbors = 5, cluster_data = FALSE, leiden_resolution = 0.5, umap_min_dist = 0.25,
                                                  report_values_as = "normalized counts", feature_reduction_method = "pca")
ggsave(filename = "flu_umap_cluster_1p5_by_type_footprint.png", plot = annotation_type_1p5, device = "png", path = "J:/U54_grant/sc/out_figures",
       width = 24, height = 22, units = "in", dpi = 600, limitsize = F, bg = "white")


# testing seutools::seurat_feature_violin_test()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")

sfvt_1 <- seutools::seurat_feature_violin_test(seurat_object = seu,
                                               plot_features = c("TRAV1-2","CD8A","CLEC12A","CD79A"),
                                               categorical_column = "cell_type",
                                               plot_categorical_types = c("MAIT/gd","Naive_CD8","CD14_Mono","CD16_Mono","B"),  # or "all"
                                               assay = "RNA",
                                               text_expansion = 1,
                                               condition = "media",
                                               condition_cat = "condition",
                                               test_cat = "age_group")


# testing seutools::seurat_mean_count_hm()
library(ggplot2)
library(Matrix)
library(seutools)
library(Seurat)
library(SeuratObject)

seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")

ref <- list('T cell' = 'CM_CD4',
            'Myeloid' = 'Mono',
            'T cell' = 'EM_CD8',
            'T cell' = 'MAIT/gd')
annocol <- c('T cell' = 'red',
             'Myeloid' = 'blue',
             'T cell' = 'red',
             'T cell' = 'red',
             'other' = 'grey')

smch <- seutools::seurat_mean_count_hm(seurat_object = seu,
                                       assay = 'RNA',
                                       cluster_column = 'cell_type',
                                       plot_clusters = c('CM_CD4','Mono','EM_CD8','MAIT/gd'),
                                       pid_column = 'pid',
                                       gene_set = c('CD8A','TNF','IL7R','IL32','SELL','LEF1','CD14','CLEC12A','TRAV1-2'),
                                       low_mid_high_cols = c("#1976D2","black","#FF1D23"),
                                       scale_per_gene = TRUE,
                                       cluster_rows = FALSE,
                                       auto_order_genes = TRUE,
                                       get_legend = FALSE,
                                       split_by_pid = TRUE,
                                       text_expansion_factor = 0.5,
                                       cluster_annotation_color = annocol, # or NULL
                                       cluster_annotation_ref = ref)
ggsave(filename = "flu_mean_hm_1.pdf", plot = smch$out_figure, device = "pdf", path = "J:/U54_grant/sc/out_figures",
       width = 5, height = 2, units = "in", dpi = 300, limitsize = F, bg = "white")


# testing seutools::seurat_dge()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")


# testing seutools::seurat_feature_violin_test()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")


cc <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/GRV_p_cluster_less4_more3/cellchat_p_GRV.rds")
cc_list <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/GRV_p_cluster_less4_more3/cellchat_p_GRV_list.rds")
lr.union <- union(cc_list[[1]]@net$LRs, cc_list[[2]]@net$LRs)
cnsh <- cellchat_netAnalysis_signalingRole_heatmap(object = cc_list[[1]], signaling = lr.union, pattern = "outgoing",
                                                   slot.name = "net", color.use = NULL, color.heatmap = "BuGn", title = NULL,
                                                   font.size.expansion = 1, cluster.rows = FALSE, cluster.cols = FALSE)


# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/GRV_p_cluster_less4_more3/cellchat_p_GRV.rds"); names(obj@var.features)
# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/GRV_pc_cluster_contr9_progr7/cellchat_pc_GRV.rds"); names(obj@var.features)
# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/Media_cl_cluster_controller11_less5/cellchat_cl_Media.rds"); names(obj@var.features)
# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/Media_p_cluster_less5_more3/cellchat_p_Media.rds"); names(obj@var.features)
# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/Media_pc_cluster_contr11_progr8/cellchat_pc_Media.rds"); names(obj@var.features)
# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/MTB300_p_cluster_less5_more3/cellchat_p_MTB300.rds"); names(obj@var.features)
# obj <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/MTB300_pc_cluster_contr11_progr8/cellchat_pc_MTB300.rds"); names(obj@var.features)

what_net <- "net"
merged_hm <- seutools::cellchat_netAnalysis_signalingRole_merged_heatmap(cellchat_object = obj, slot.name = what_net, color.use = NULL,
                                                                         name_vector_key = c("progr_less_4"="progressor <12mo","progr_more_3"="progressor",
                                                                                             "controller_9"="controller","progressor_7"="progressor",
                                                                                             "controller_11"="controller","progr_less_5"="progressor <12mo",
                                                                                             "progressor_8"="progressor"),
                                                                         font.size.expansion = ifelse(what_net=="net",0.5,0.8))
pdf(file = paste0("J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/MTB300_pc_cluster_contr11_progr8/cellchat_pc_MTB300_signaling_merge_",what_net,".pdf"),
    width = ifelse(what_net=="net",16,16), height = ifelse(what_net=="net",15,10))
merged_hm
dev.off()
