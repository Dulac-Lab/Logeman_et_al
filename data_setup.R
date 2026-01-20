library(Seurat)
library(Matrix)
library(ArchR)
library(tidyverse)

addArchRThreads(threads = 12)
addArchRGenome("mm10")
setwd("/Volumes/seqSSD/01_Logeman")

load_samples_RNA <- function(mex_dir, feature_column = 2, prefix_cells = TRUE, set_idents_from = "Cluster") {
  dirs <- list.dirs(mex_dir, recursive = FALSE, full.names = TRUE)
  
  xs <- lapply(dirs, function(d) {
    s <- basename(d); fp <- function(x) file.path(d, paste0(s, "_", x))
    m  <- as(readMM(fp("matrix.mtx")), "CsparseMatrix")
    f  <- read.delim(fp("features.tsv"), header = FALSE, sep = "\t", check.names = FALSE)[[feature_column]]
    bc <- readLines(fp("barcodes.tsv"))
    md <- read.delim(fp("metadata.tsv"), sep = "\t", check.names = FALSE)
    
    rownames(m) <- make.unique(f); colnames(m) <- bc
    rownames(md) <- md$barcode; md$barcode <- NULL
    
    if (prefix_cells) {
      cn <- paste0(s, "_", colnames(m))
      colnames(m) <- cn
      rownames(md) <- cn
    }
    
    md <- md[colnames(m), , drop = FALSE]; md$orig.ident <- s
    list(m = m, md = md)
  })
  
  m_all  <- Reduce(cbind2, lapply(xs, `[[`, "m"))
  md_all <- do.call(rbind, lapply(xs, `[[`, "md"))
  
  obj <- CreateSeuratObject(m_all)
  obj <- AddMetaData(obj, md_all[colnames(obj), , drop = FALSE])
  if (!is.null(set_idents_from) && set_idents_from %in% colnames(obj[[]]))
    Idents(obj) <- factor(obj[[set_idents_from]][, 1])
  obj
}

combined_RNA <- load_samples_RNA("RNA/input_data/RNA_count_matrix")
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
  #outputNames  = outputNames,
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
  #geneAnnotation = genes,
  #genomeAnnotation = genome,
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

