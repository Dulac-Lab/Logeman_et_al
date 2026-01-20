library(ArchR)
library(Seurat)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

#combined_ArchR <- readRDS("/Volumes/seqSSD/combined_ArchR/combined_ArchR_01_naive/Save-ArchR-Project.rds")

setwd("/Volumes/seqSSD/ATAC")
md <- read.csv("/Users/brandonlogeman/Desktop/ATAC/cell_metadata.csv", row.names = 1, header = TRUE)

# Define gene ranges
genes <- createGeneAnnotation(
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  OrgDb = org.Mm.eg.db
)

# Define genome
genome <- createGenomeAnnotation(
  genome = BSgenome.Mmusculus.UCSC.mm10
)

Gal_Fathers_1 <- "ArrowFiles/Gal_Fathers_1.arrow"
Gal_Mothers_1 <- "ArrowFiles/Gal_Mothers_1.arrow"
Gal_VirginFemales_1 <- "ArrowFiles/Gal_VirginFemales_1.arrow"
Gal_VirginMales_1 <- "ArrowFiles/Gal_VirginMales_1.arrow"

Brs3_Fathers_1 <- "ArrowFiles/Brs3_Fathers_1.arrow"
Brs3_Fathers_2 <- "ArrowFiles/Brs3_Fathers_2.arrow"
Brs3_Mothers_1 <- "ArrowFiles/Brs3_Mothers_1.arrow"
Brs3_Mothers_2 <- "ArrowFiles/Brs3_Mothers_2.arrow"
Brs3_VirginFemales_1 <- "ArrowFiles/Brs3_VirginFemales_1.arrow"
Brs3_VirginFemales_2 <- "ArrowFiles/Brs3_VirginFemales_2.arrow"
Brs3_VirginMales_1 <- "ArrowFiles/Brs3_VirginMales_1.arrow"
Brs3_VirginMales_2 <- "ArrowFiles/Brs3_VirginMales_2.arrow"

Ucn3_Fathers_1 <- "ArrowFiles/Ucn3_Fathers_1.arrow"
Ucn3_Fathers_2 <- "ArrowFiles/Ucn3_Fathers_2.arrow"
Ucn3_Fathers_3 <- "ArrowFiles/Ucn3_Fathers_3.arrow"
Ucn3_Mothers_1 <- "ArrowFiles/Ucn3_Mothers_1.arrow"
Ucn3_Mothers_2 <- "ArrowFiles/Ucn3_Mothers_2.arrow"
Ucn3_Mothers_3 <- "ArrowFiles/Ucn3_Mothers_3.arrow"
Ucn3_VirginFemales_1 <- "ArrowFiles/Ucn3_VirginFemales_1.arrow"
Ucn3_VirginFemales_2 <- "ArrowFiles/Ucn3_VirginFemales_2.arrow"
Ucn3_VirginFemales_3 <- "ArrowFiles/Ucn3_VirginFemales_3.arrow"
Ucn3_VirginMales_1 <- "ArrowFiles/Ucn3_VirginMales_1.arrow"
Ucn3_VirginMales_2 <- "ArrowFiles/Ucn3_VirginMales_2.arrow"
Ucn3_VirginMales_3 <- "ArrowFiles/Ucn3_VirginMales_3.arrow"


# Merge arrow files into single project
combined_ArchR <- ArchRProject(
  ArrowFiles = c(Gal_Fathers_1, Gal_Mothers_1, Gal_VirginFemales_1, Gal_VirginMales_1, Brs3_Fathers_1, Brs3_Fathers_2, Brs3_Mothers_1, Brs3_Mothers_2, Brs3_VirginFemales_1, Brs3_VirginFemales_2, Brs3_VirginMales_1, Brs3_VirginMales_2, Ucn3_Fathers_1, Ucn3_Fathers_2, Ucn3_Fathers_3, Ucn3_Mothers_1, Ucn3_Mothers_2, Ucn3_Mothers_3, Ucn3_VirginFemales_1, Ucn3_VirginFemales_2, Ucn3_VirginFemales_3, Ucn3_VirginMales_1, Ucn3_VirginMales_2, Ucn3_VirginMales_3),
  outputDirectory = 'TEST',
  copyArrows = FALSE,
  geneAnnotation = genes,
  genomeAnnotation = genome,
  showLogo = FALSE
)

rm(Gal_Fathers_1, Gal_Mothers_1, Gal_VirginFemales_1, Gal_VirginMales_1, Brs3_Fathers_1, Brs3_Fathers_2, Brs3_Mothers_1, Brs3_Mothers_2, Brs3_VirginFemales_1, Brs3_VirginFemales_2, Brs3_VirginMales_1, Brs3_VirginMales_2, Ucn3_Fathers_1, Ucn3_Fathers_2, Ucn3_Fathers_3, Ucn3_Mothers_1, Ucn3_Mothers_2, Ucn3_Mothers_3, Ucn3_VirginFemales_1, Ucn3_VirginFemales_2, Ucn3_VirginFemales_3, Ucn3_VirginMales_1, Ucn3_VirginMales_2, Ucn3_VirginMales_3)
rm(genes, genome)

# add each column of metadata
common <- intersect(rownames(md), combined_ArchR$cellNames)
for (nm in colnames(md)) {
  combined_ArchR <- addCellColData(
    ArchRProj = combined_ArchR,
    data  = md[[nm]],
    cells = common,
    name  = nm,
    force = TRUE   # overwrite if column exists :contentReference[oaicite:1]{index=1}
  )
}
rm(md, common, nm)

setwd("/Volumes/seqSSD/ATAC")
clusters_to_make <- c("I-14", "E-1")

for (cl in clusters_to_make) {
  
  # folder for this cluster
  out_dir <- file.path(getwd(), cl)
  
  # create folder if needed
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("Directory created: ", out_dir)
  } else {
    message("Directory already exists: ", out_dir)
  }
  
  # cells for this cluster
  cells_keep <- getCellNames(combined_ArchR)[combined_ArchR$Cluster %in% cl]
  
  # skip if no cells found
  if (length(cells_keep) == 0) {
    warning("No cells found for cluster: ", cl, " (skipping)")
    next
  }
  
  # subset project
  sub_proj <- subsetArchRProject(
    ArchRProj = combined_ArchR,
    cells = cells_keep,
    outputDirectory = out_dir,
    force = TRUE
  )
  
  message("Done cluster ", cl, " | cells: ", length(cells_keep))
}
rm(clusters_to_make, cells_keep, cl, clusters_to_make, out_dir)

# Perform iterative Latent Semantic Indexing
combined_ArchR <- addIterativeLSI(
  ArchRProj = combined_ArchR,
  useMatrix = "TileMatrix", # either PeakMatrix or TileMatrix
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    maxClusters = 6,
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  UMAPParams = list(n_neighbors = 40,
                    min_dist = 0.4,
                    metric = "cosine",
                    verbose = FALSE,
                    fast_sgd = TRUE
  ),
)

#Perform UMAP
combined_ArchR <- addUMAP(
  ArchRProj = combined_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

p1 <- plotEmbedding(
  ArchRProj = combined_ArchR,
  colorBy = "cellColData",
  name = "Cluster",
  embedding = "UMAP")
p1