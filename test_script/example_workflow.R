# example workflow
library(seutools)

und_cl <- paste0("Undecided",1:5)

# testing seutools::seurat_mast()
und_seu_mast <- seurat_mast(seurat_object = seu, freq_expressed = 0.1, fc_threshold = log2(1.5),
                            test_per_cluster = TRUE, test_clusters = und_cl, cluster_column = "cell_type",
                            test_per_category = FALSE, test_categories = c("younger","older"),
                            category_column = "age_group", test_per_condition = FALSE,
                            test_condition = "all", condition_column = "condition", pid_column = "pid")


# testing seutools::seurat_feature_overlay()
library(Matrix)
library(seutools)
library(Seurat)

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
