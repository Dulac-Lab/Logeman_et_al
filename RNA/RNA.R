library(Seurat)
library(tidyverse)
library(Augur)

setwd("/Volumes/seqSSD/01_Logeman")

IEG.genes = c("Fos","Npas4","Arc","Jun","Egr1",'Junb','Fosb','Nr4a1', "Fosl2", "Nr4a3", "Egr4")
Sex_linked_genes <- c("Xist", "Tsix", "Eif2s3y", "Uty", "Ddx3y", "Kdm5d", "Hya", "Gm21860", "Gm47283", "Gm29650", "Usp9y")
comparisons = list(
  c("VirginMale", "VirginMale"),
  c("Father", "Father"),
  c("VirginFemale", "VirginFemale"),
  c("Mother", "Mother"),
  c("VirginMale", "Father"),
  c("VirginMale", "Mother"),
  c("VirginMale", "VirginFemale"),
  c("Father", "Mother"),
  c("Father", "VirginFemale"),
  c("Mother", "VirginFemale")
)

merfish <- readRDS(file.path(getwd(), "RNA/Results/all_clusters/MERFISH.rds"))

merfish <- SCTransform(
  object = merfish,
  assay = "RNA",
  conserve.memory = TRUE
)

RNA <- readRDS(file.path(getwd(), "RNA/Results/all_clusters/all_clusters.rds"))

RNA <- SCTransform(
  object = RNA,
  assay = "RNA",
  conserve.memory = TRUE
)

sample.features <- setdiff(VariableFeatures(object = RNA), c(IEG.genes, Sex_linked_genes))

RNA <- RunPCA(
  object = RNA,
  assay = "SCT",
  features = sample.features,
  npcs = 150,
  verbose = TRUE
)

RNA <- RunUMAP(
  object = RNA,
  dims = 1:150,
  reduction = "pca"
)

rm(sample.features)

Integrated.anchors <- FindTransferAnchors(reference = merfish,
                                          query = RNA,
                                          normalization.method = "SCT",
                                          query.assay = "SCT",
                                          reduction = "cca")

Integrated.predictions <- TransferData(anchorset = Integrated.anchors,
                                       refdata = merfish$Neuron_cluster_ID,
                                       weight.reduction = "cca")

predictions <- Integrated.predictions %>%
  rownames_to_column('barcode') %>%
  select('barcode', contains('prediction.score')) %>%
  select(-prediction.score.max) %>%
  full_join(y=FetchData(
    object = RNA,
    vars = c('Cluster')
  ) %>% rownames_to_column('barcode'),
  by = 'barcode') %>%
  column_to_rownames('barcode') %>%
  group_by(cluster) %>%
  summarise_all(mean) %>%
  column_to_rownames('Cluster') %>%
  rename_with(~str_remove(., 'prediction.score.')) %>%
  select(order(colnames(.)))

rm(MERFISH, Integrated.anchors, Integrated.predictions)


### Train random forest classifier to predict State for each cluster
DefaultAssay(object = RNA) <- "RNA"
RNA <- NormalizeData(object = RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA <- ScaleData(object = RNA)

############################# Calculate ###########################################
classifier.results.comparision.list <- list()
for (seurat_cluster in unique(RNA$Cluster)) {
  list.name <- paste('classifier.results', seurat_cluster, sep = '.cluster_')
  Idents(RNA) <- "Cluster"
  cluster.subset <- subset(RNA, idents = seurat_cluster)
  # Create a list to store Augur results
  classifier_results = list()
  # Run Augur on each comparison
  for (i in seq_along(comparisons)) {
    Idents(cluster.subset) = "State"
    if (length(unique(comparisons[[i]])) > 1){
      # Subset the Seurat object to these two experimental conditions only
      sub <- subset(cluster.subset, idents = comparisons[[i]])
    } else {
      # Duplicate Seurat object 
      sub1 <- subset(cluster.subset, idents = unique(comparisons[[i]]))
      sub1$State <- paste(sub1$State, "1", sep = "_")
      sub2 <- subset(cluster.subset, idents = unique(comparisons[[i]]))
      sub2$State <- paste(sub2$State, "2", sep = "_")
      sub <- merge(sub1,
                   y = sub2,
                   add.cell.ids = c('dup_1', 'dup_2'),
                   merge.data = TRUE
      )
    }
    # Subset the most variable features minus IEGs and Sex linked genes
    gene.list <- rownames(cluster.subset)
    gene.list <- setdiff(gene.list, c(IEG.genes, Sex_linked_genes))
    sub <- subset(sub, features = gene.list)
    sub <- FindVariableFeatures(sub, selection.method = 'vst', nfeatures = 100)
    var.genes <- sub@assays[["RNA"]]@var.features
    expr <- GetAssayData(sub)[var.genes,]
    meta <- sub@meta.data
    
    # Run Augur
    augur = calculate_auc(
      expr,
      meta = meta,
      label_col = 'State',
      cell_type_col = 'Cluster',
      var_quantile = 1,
      feature_perc = 1,
      rf_params = list(
        trees = 100,
        mtry = 2,
        min_n = NULL,
        importance = "accuracy")
    )
    
    # Store the results
    comparison_name = paste(comparisons[[i]], collapse = ":")
    classifier_results[[comparison_name]] = augur
  }
  classifier.results.comparision.list[[list.name]] <- classifier_results
}
rm(augur, expr, meta, sub, sub1, sub2, var.genes, i, comparison_name, gene.list, classifier_results)

# Save results
saveRDS(classifier.results.comparision.list, file = 'RNA/Results/all_clusters/classifier_results.rds')

