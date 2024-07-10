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


# testing seutools::seurat_dge()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")


# testing seutools::seurat_feature_violin_test()
seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")


#################
