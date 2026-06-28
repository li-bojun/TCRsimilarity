library(SoupX)
library(Matrix)
library(Seurat)

# from R console: system("Rscript runsoupx.r your_data_dirs.csv")

# csv expected format:
# raw, filtered, output
# /path/raw1, /path/filtered1, /path/output1
# /path/raw2, /path/filtered2, /path/output2

run_soupx <- function(input_dir_raw, input_dir_filtered, output_dir){
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  raw <- Read10X(input_dir_raw)
  filtered <- Read10X(input_dir_filtered)
  
  if (is.list(raw)) raw <- raw[["Gene Expression"]]
  if (is.list(filtered)) filtered <- filtered[["Gene Expression"]]
  
  sc <- SoupChannel(tod = raw, toc = filtered)
  
  seu <- CreateSeuratObject(filtered)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu)
  
  clusters <- Idents(seu)
  sc <- setClusters(sc, clusters)
  
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  
  writeMM(out, file.path(output_dir, "matrix.mtx"))
  
  write.table(colnames(out),
              gzfile(file.path(output_dir, "barcodes.tsv.gz")),
              quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = FALSE)
  
  write.table(
    data.frame(
      gene_id = rownames(out),
      gene_name = rownames(out),
      feature_type = "Gene Expression"
    ),
    gzfile(file.path(output_dir, "features.tsv.gz")),
    quote = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE
  )
  
  system(paste("gzip -f", file.path(output_dir, "matrix.mtx")))
  
  cat("Sample finished: ", output_dir, "\n")
}

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1){
  stop("Please provide the path to the samples CSV as the first argument.")
}

samples_csv <- args[1]

samples <- read.csv(samples_csv, stringsAsFactors = FALSE)
for(i in 1:nrow(samples)){
  cat("Processing sample", i, "of", nrow(samples), "\n")
  run_soupx(samples$raw[i],
            samples$filtered[i],
            samples$output[i])
}
print("All SoupX runs complete.")
