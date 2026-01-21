library(Seurat)
library(Matrix)
library(ArchR)
library(tidyverse)

addArchRThreads(threads = 10)
addArchRGenome("mm10")
setwd("/Volumes/seqSSD/01_Logeman")

load_samples_RNA <- function(mex_dir, feature_column = 2, prefix_cells = TRUE, set_idents_from = "Cluster") {
  d <- normalizePath(mex_dir, winslash = "/", mustWork = TRUE)
  pfx <- sort(sub("_barcodes\\.tsv$", "", list.files(d, "_barcodes\\.tsv$", full.names = FALSE)))
  stopifnot(length(pfx) > 0)
  
  xs <- lapply(pfx, function(s) {
    fp <- function(x) file.path(d, paste0(s, "_", x))
    m  <- methods::as(Matrix::readMM(fp("matrix.mtx")), "CsparseMatrix")
    rn <- make.unique(read.delim(fp("features.tsv"), header = FALSE, sep = "\t", check.names = FALSE)[[feature_column]])
    bc <- readLines(fp("barcodes.tsv"))
    md <- read.delim(fp("metadata.tsv"), sep = "\t", check.names = FALSE)
    
    stopifnot(nrow(m) == length(rn), ncol(m) == length(bc))
    rownames(m) <- rn; colnames(m) <- bc
    rownames(md) <- md$barcode; md$barcode <- NULL; md <- md[bc, , drop = FALSE]
    if (prefix_cells) { cn <- paste0(s, "_", bc); colnames(m) <- cn; rownames(md) <- cn }
    md$orig.ident <- s
    list(m = m, md = md)
  })
  
  genes <- unique(unlist(lapply(xs, function(z) rownames(z$m)), use.names = FALSE))
  
  reidx <- function(m) {
    map <- match(rownames(m), genes)
    Matrix::sparseMatrix(
      i = map[m@i + 1L], j = rep(seq_len(ncol(m)), diff(m@p)), x = m@x,
      dims = c(length(genes), ncol(m)), dimnames = list(genes, colnames(m))
    )
  }
  
  m_all <- Reduce(Matrix::cbind2, lapply(xs, function(z) reidx(z$m)))
  md_all <- do.call(rbind, lapply(xs, `[[`, "md"))
  
  obj <- Seurat::CreateSeuratObject(m_all)
  obj <- Seurat::AddMetaData(obj, md_all[colnames(obj), , drop = FALSE])
  if (!is.null(set_idents_from) && set_idents_from %in% colnames(obj[[]]))
    Seurat::Idents(obj) <- factor(obj[[set_idents_from]][, 1])
  obj
}

combined_RNA <- load_samples_RNA("RNA/input_data/RNA_count_matrix/")
dir.create("RNA/Results/all_clusters", recursive = TRUE, showWarnings = FALSE)
saveRDS(combined_RNA, "RNA/Results/all_clusters/all_clusters.rds")

### SUB CLUSTER ###

## Seurat version of “split by cluster name + save each subset in its own folder”
clusters_to_make <- c("I-14", "E-1")
for (cl in clusters_to_make) {
  
  # folder for this cluster
  out_dir <- file.path("RNA/Results", cl)
  
  # create folder if needed
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("Directory created: ", out_dir)
  } else {
    message("Directory already exists: ", out_dir)
  }
  
  # subset object
  Idents(combined_RNA) <- "Cluster"
  sub_seu <- subset(combined_RNA, idents = cl)
  
  # save (RDS)
  safe_cl <- gsub("[^A-Za-z0-9._-]", "_", cl)
  saveRDS(sub_seu, file = file.path(out_dir, paste0(safe_cl, ".rds")))
  
  message("Done cluster ", cl)
}

rm(clusters_to_make, cl, out_dir, sub_seu, safe_cl)

#####################################################
#####################################################
#####################################################

# Import data
merfish_dat <- read_csv(file.path(getwd(), "RNA/input_data/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv"))

# Assign physiological state to samples
merfish_dat <- merfish_dat %>%
  mutate(State = case_when(
    Animal_sex == "Male" & Behavior == "Parenting" ~ "Father",
    Animal_sex == "Male" & Behavior == "Aggression to pup" ~ "VirginMale",
    Animal_sex == "Male" & Behavior == "Aggression to adult" ~ "VirginMale",
    Animal_sex == "Male" & Behavior == "Mating" ~ "VirginMale",
    Animal_sex == "Male" & Behavior == "Naive" ~ "VirginMale",
    Animal_sex == "Female" & Behavior == "Parenting" ~ "Mother",
    Animal_sex == "Female" & Behavior == "Virgin Parenting" ~ "VirginFemale",
    Animal_sex == "Female" & Behavior == "Mating" ~ "VirginFemale",
    Animal_sex == "Female" & Behavior == "Naive" ~ "VirginFemale"))

merfish_dat <- merfish_dat %>%
  dplyr::select(1:3, State, everything(), -Blank_1, -Blank_2, -Blank_3, -Blank_4, -Blank_5, -Fos)

#Build counts matrix: genes (rows) x cells (cols)
counts_mat <- merfish_dat %>%
  select(Cell_ID, Ace2:Vgf) %>%
  column_to_rownames("Cell_ID") %>%
  as.matrix()

# Convert to sparse and transpose: (cells x genes) -> (genes x cells)
counts_mat <- t(Matrix(counts_mat, sparse = TRUE))

#Build metadata once: cells (rows) x metadata (cols)
meta_df <- merfish_dat %>%
  select(
    Cell_ID, Neuron_cluster_ID, Animal_ID, Animal_sex, Behavior,
    Cell_class, State, Bregma, Centroid_X, Centroid_Y
  ) %>%
  distinct(Cell_ID, .keep_all = TRUE) %>%
  column_to_rownames("Cell_ID")

#Create Seurat object with metadata attached
MERFISH <- CreateSeuratObject(
  counts    = counts_mat,
  project   = "MERFISH",
  assay     = "RNA",
  meta.data = meta_df
)

# Select relevant clusters
MERFISH <- subset(MERFISH, subset = Neuron_cluster_ID %in% unique(combined_RNA@meta.data$Cluster))

# Save RDS
saveRDS(MERFISH, "RNA/Results/all_clusters/MERFISH.rds")

#####################################################
#####################################################
#####################################################

# Find all .tsv.gz fragment files
fragPath <- normalizePath(file.path(getwd(), "ATAC/input_data/FragmentFiles"), mustWork = FALSE)
inputFiles <- list.files(
  path = fragPath,
  pattern = "\\.tsv\\.gz$",
  full.names = TRUE
)

# Sample names = file name
sampleNames <- sub("\\.tsv\\.gz$", "", basename(inputFiles))

# Where to write Arrow files
outDir <- file.path(getwd(), "ATAC", "Results", "all_clusters", "ArrowFiles")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# Create Arrow files
setwd(outDir)
arrowFiles <- createArrowFiles(
  inputFiles   = inputFiles,
  sampleNames  = sampleNames,
  QCDir = file.path(getwd(), "QualityControl"),
  minTSS = 1,
  minFrags = 10,
  maxFrags = 1000000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
normalizePath(file.path(getwd(), arrowFiles), mustWork = FALSE)
setwd("../../../..")

### DELETE ME ###

arrowFiles <- list.files(
  path = "/Volumes/seqSSD/01_Logeman/ATAC/Results/all_clusters/ArrowFiles",
  pattern = ".arrow$",
  full.names = TRUE
)

#################

# Merge arrow files into single project
combined_ArchR <- ArchRProject(
  ArrowFiles = arrowFiles,
  outputDirectory = 'ATAC/Results/all_clusters',
  copyArrows = FALSE,
  showLogo = FALSE
)

# add each column of metadata
md <- do.call(rbind, lapply(
  list.files(fragPath, pattern="_metadata\\.csv$", full.names=TRUE),
  read.csv, stringsAsFactors=FALSE
))
rownames(md) <- md$barcode
md$barcode <- NULL

common <- intersect(rownames(md), combined_ArchR$cellNames)
for (nm in colnames(md)) {
  combined_ArchR <- addCellColData(
    ArchRProj = combined_ArchR,
    data  = md[[nm]],
    cells = common,
    name  = nm,
    force = TRUE
  )
}
rm(md, common, nm, arrowFiles, fragPath, inputFiles, outDir, sampleNames)

combined_ArchR <- saveArchRProject(
  ArchRProj = combined_ArchR,
  outputDirectory = 'ATAC/Results/all_clusters',
  overwrite = TRUE,
  load = TRUE
  )

### SUB CLUSTER ###

clusters_to_make <- c("I-14", "E-1")
for (cl in clusters_to_make) {
  
  # folder for this cluster
  out_dir <- file.path("ATAC/Results", cl)
  
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
  
  saveArchRProject(
    ArchRProj = sub_proj,
    outputDirectory = out_dir,
    overwrite = TRUE,
    load = FALSE
  )
  
  message("Done cluster ", cl, " | cells: ", length(cells_keep))
}
rm(clusters_to_make, cells_keep, cl, clusters_to_make, out_dir)

