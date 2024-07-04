
library(Seurat)
library(Matrix)

if(F) {
  seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object.rds")
  page <- read.csv(file = "J:/U54_grant/sc/inputs/pid_age.csv", check.names = FALSE)
  if(mean(unique(seu@meta.data$pid) %in% page$pid)!=1) {
    stop("one or more pid not found in age metadata")
  } else {
    ages <- page$age; names(ages) <- page$pid
    seu <- AddMetaData(object = seu, metadata = as.numeric(ages[seu@meta.data$pid]), col.name = "age")
    seu <- AddMetaData(object = seu, metadata = ifelse(seu$age>55, "older", "younger"), col.name = "age_group")
  }
}

if(F) {
  adt_ct <- read.csv(file = "J:/U54_grant/sc/inputs/adt_adata_counts.csv", check.names = FALSE, header = FALSE)
  adt_obs <- read.csv(file = "J:/U54_grant/sc/inputs/adt_adata_obs.csv", check.names = FALSE)
  adt_var <- read.csv(file = "J:/U54_grant/sc/inputs/adt_adata_var.csv", check.names = FALSE)

  adt_ct <- t(adt_ct)
  row.names(adt_ct) <- row.names(adt_var)
  colnames(adt_ct) <- adt_obs$barcode_1

  seu_adt <- Seurat::CreateSeuratObject(counts = adt_ct, assay = "ADT", meta.data = adt_obs)
  seu_adt@assays[["ADT"]]@layers[["data"]] <- seu_adt@assays[["ADT"]]@layers[["counts"]]

  seu_adt <- AddMetaData(object = seu_adt, metadata = as.numeric(ages[seu_adt@meta.data$pid]), col.name = "age")
  seu_adt <- AddMetaData(object = seu_adt, metadata = ifelse(seu_adt$age>55, "older", "younger"), col.name = "age_group")

  adata_map <- read.csv(file = "J:/U54_grant/sc/inputs/adata_umap.csv", check.names = FALSE)
  adata_obs <- read.csv(file = "J:/U54_grant/sc/inputs/adata_obs.csv", check.names = FALSE)
  map1 <- adata_map$UMAP1; names(map1) <- gsub("-1_","_",adata_obs$barcode)
  map2 <- adata_map$UMAP2; names(map2) <- gsub("-1_","_",adata_obs$barcode)

  adt_umap <- data.frame(UMAP1 = map1[seu_adt@meta.data$barcode_1], UMAP2 = map2[seu_adt@meta.data$barcode_1]); adt_umap <- as.matrix(adt_umap)
  row.names(adt_umap) <- seu_adt@meta.data$barcode_1
  adt_dr <- Seurat::CreateDimReducObject(embeddings = adt_umap, key = "UMAP_", assay = "ADT")

  seu_adt[['umap']] <- adt_dr

  saveRDS(object = seu_adt, file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")
} else if(F){
  seu_adt <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_adt_seurat_object.rds")
}

if(F) {
  # und_ct <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_counts.csv", check.names = FALSE, header = FALSE)
  und_ct <- Matrix::readMM(file = "J:/U54_grant/sc/inputs/undecided_adata_counts.mtx")
  und_obs <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_obs.csv", check.names = FALSE)
  und_var <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_var.csv", check.names = FALSE)

  und_latent <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_latent.csv")
  und_umap <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_umap.csv")
  und_norm <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_scvi_norm_counts.csv")

  und_ct <- Matrix::t(und_ct)
  row.names(und_ct) <- und_var[,1]
  colnames(und_ct) <- und_obs$barcode

  seu_und <- Seurat::CreateSeuratObject(counts = und_ct, assay = "RNA", meta.data = und_obs)

  seu_und <- AddMetaData(object = seu_und, metadata = as.numeric(ages[seu_und@meta.data$pid]), col.name = "age")
  seu_und <- AddMetaData(object = seu_und, metadata = ifelse(seu_und$age>55, "older", "younger"), col.name = "age_group")

  adata_map <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_umap.csv", check.names = FALSE)
  adata_obs <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_obs.csv", check.names = FALSE)
  map1 <- adata_map$UMAP1; names(map1) <- gsub("-1_","_",adata_obs$barcode)
  map2 <- adata_map$UMAP2; names(map2) <- gsub("-1_","_",adata_obs$barcode)

  und_umap <- data.frame(UMAP1 = map1[seu_und@meta.data$barcode], UMAP2 = map2[seu_und@meta.data$barcode]); und_umap <- as.matrix(und_umap)
  row.names(und_umap) <- seu_und@meta.data$barcode
  und_dr <- Seurat::CreateDimReducObject(embeddings = und_umap, key = "UMAP_", assay = "RNA")

  seu_und[['umap']] <- und_dr

  saveRDS(object = seu_und, file = "J:/U54_grant/sc/inputs/pbmc_flu_undecided_seurat_object.rds")
} else if(F){
  seu_und <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_undecided_seurat_object.rds")
}

if(F) {
  seu_all <- merge(x = seu, y = seu_und)
  seu_all <- JoinLayers(object = seu_all)
  seu_all_ct <- Seurat::GetAssayData(object = seu_all, assay = "RNA", layer = "counts")
  seu_all_obj <- Seurat::CreateSeuratObject(counts = seu_all_ct, assay = "RNA", meta.data = seu_all@meta.data)
  seu_all_obj <- Seurat::NormalizeData(object = seu_all_obj)

  adata_map <- read.csv(file = "J:/U54_grant/sc/inputs/adata_umap.csv", check.names = FALSE)
  adata_obs <- read.csv(file = "J:/U54_grant/sc/inputs/adata_obs.csv", check.names = FALSE)
  map1_1 <- adata_map$UMAP1; names(map1_1) <- gsub("-1_","_",adata_obs$barcode)
  map2_1 <- adata_map$UMAP2; names(map2_1) <- gsub("-1_","_",adata_obs$barcode)

  adata_map <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_umap.csv", check.names = FALSE)
  adata_obs <- read.csv(file = "J:/U54_grant/sc/inputs/undecided_adata_obs.csv", check.names = FALSE)
  map1_2 <- adata_map$UMAP1; names(map1_2) <- gsub("-1_","_",adata_obs$barcode)
  map2_2 <- adata_map$UMAP2; names(map2_2) <- gsub("-1_","_",adata_obs$barcode)

  map1 <- c(map1_1, map1_2); map2 <- c(map2_1, map2_2)

  umap_mat <- as.matrix(data.frame(UMAP1 = map1[gsub("-1_","_",seu_all_obj$barcode)], UMAP2 = map2[gsub("-1_","_",seu_all_obj$barcode)]))
  row.names(umap_mat) <- gsub("_","-1_",row.names(umap_mat))

  all_dr <- Seurat::CreateDimReducObject(embeddings = umap_mat, key = "UMAP_", assay = "RNA")
  seu_all_obj[['umap']] <- all_dr

  saveRDS(object = seu_all_obj, file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
} else if(F){
  seu <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_seurat_object_all.rds")
  seu_small <- subset(x = seu, subset = lane %in% c("9569-JO-1_multi","9569-JO-2_multi"))
  saveRDS(object = seu_small, file = "J:/U54_grant/sc/inputs/pbmc_flu_small_seurat_object_all.rds")
} else {
  library(Seurat)
  library(Matrix)
  seu_small <- readRDS(file = "J:/U54_grant/sc/inputs/pbmc_flu_small_seurat_object_all.rds")
}

if(F) {
  seu_act_nk <- subset(x = seu, subset = cell_type == "Activ NK")
  seu_mono <- subset(x = seu, subset = cell_type == "Mono")
  seu_test_dge <- merge(x = seu_act_nk, y = seu_mono)
  seu_test_dge <- JoinLayers(object = seu_test_dge, assay = "RNA", layers = "data")
  seu_test_dge <- JoinLayers(object = seu_test_dge, assay = "RNA", layers = "counts")
}

# liana_in
# pseudobulk_in
