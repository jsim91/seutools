# example workflow
library(seutools)

und_cl <- paste0("Undecided",1:5)

# goal: test Undecided clusters against other without subsetting by condition
und_seu_mast <- seurat_mast(seurat_object = seu, freq_expressed = 0.1, fc_threshold = log2(1.5),
                            test_per_cluster = TRUE, test_clusters = und_cl, cluster_column = "cell_type",
                            test_per_category = FALSE, test_categories = c("younger","older"),
                            category_column = "age_group", test_per_condition = FALSE,
                            test_condition = "all", condition_column = "condition", pid_column = "pid")
