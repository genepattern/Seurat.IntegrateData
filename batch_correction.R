suppressMessages(suppressWarnings(library("log4r")))
logfile <- "batch_correction_log.txt"
console_appender <- console_appender(layout = default_log_layout())
file_appender <- file_appender(logfile, append = TRUE, layout = default_log_layout())
logger <- log4r::logger(threshold = "INFO", appenders = list(console_appender, file_appender))

info(logger, message = "==========================================================")
info(logger, message = "Loading libraries: Seurat, optparse, log4r")
suppressMessages(suppressWarnings(library("Seurat")))
suppressMessages(suppressWarnings(library("optparse")))
info(logger, message = "Loaded libraries: Seurat, optparse, log4r")
info(logger, message = "==========================================================")

# # ====================================
# # PARAMETERS for PCA
# ncomps = 30 (Default)
# # ====================================
# # PARAMETERS for Violin Plot
# nCountRNA: TRUE (Default)
# nFeatureRNA: TRUE (Default)
# # ====================================
# # PARAMETERS for OUTPUT
# output_file_name: "batch_processing_results" (Default)
# # ====================================

# Parse Input Arguments
parser = OptionParser()
# ====================================
parser <- add_option(parser, c("--file_list"), help = "List of files to load. Separate file names with '\ ' ")
# ====================================
# PARAMETERS for PCA
parser <- add_option(parser, c("--ncomps"), type = "integer", default = 30, help = "How many PCA components to use (?) (default = 30)")
# ====================================
# PARAMETERS for Violin Plot
parser <- add_option(parser, c("--nCount_RNA"), type = "logical", default = TRUE, help = "Display nCount_RNA feature on Violin Plot (default = TRUE)")
parser <- add_option(parser, c("--nFeature_RNA"), type = "logical", default = TRUE, help = "Display nFeature_RNA feature on Violin Plot (default = TRUE)")
# ====================================
# PARAMETERS for Output File
parser <- add_option(parser, c("--output_file_name"), default = "batch_processing_results", help = "Output File Name for Batch Correction Analysis (default = 'batch_processing_results')")
# ====================================
info(logger, message = "==========================================================")
args <- parse_args(parser)
info(logger, message = "Parameters used:")
info(logger, message = paste("help:", args$help))
info(logger, message = paste("file_list:", args$file_list))
info(logger, message = paste("ncomps:", args$ncomps))
info(logger, message = paste("nCount_RNA:", args$nCount_RNA))
info(logger, message = paste("nFeature_RNA:", args$nFeature_RNA))
info(logger, message = paste("output_file_name:", args$output_file_name))
info(logger, message = "==========================================================")
# ====================================

con <- file(args$file_list, open = "r")
lines = readLines(con)
data_list <- list()
files <- list()
condensed_file_names <- list()

info(logger, message = "==========================================================")
info(logger, message = paste("Reading Files:", args$file_list))

for(i in 1:length(lines)){
  data_list[[i]] <- read.table(file = lines[[i]], header = TRUE, sep = "\t", row.names = 1)
  files[[i]] <- lines[[i]]
  condensed_file_names[[i]] <- tail(strsplit(lines[[i]], "/")[[1]], 1)
}

info(logger, message = paste("Read Files:", args$file_list))
info(logger, message = "==========================================================")

seurat_objects <- list()

info(logger, message = "==========================================================")
info(logger, message = paste("Creating Seurat Objects for:", args$file_list))

for (i in 1:length(data_list)){
  seurat_objects[[i]] <- CreateSeuratObject(counts = data_list[[i]], project = condensed_file_names[[i]])
  if (length(levels(seurat_objects[[i]])) > 1){
    levels(seurat_objects[[i]]@active.ident) <- rep(condensed_file_names[[i]], length(levels(seurat_objects[[i]])))
  }
}


info(logger, message = paste("Created Seurat Objects:", args$file_list))
info(logger, message = "==========================================================")

info(logger, message = "==========================================================")
info(logger, message = paste("Normalizing Data and Finding Variable Features on Seurat Objects"))
for(i in 1:length(seurat_objects)){
  seurat_objects[[i]] <- NormalizeData(seurat_objects[[i]], verbose = FALSE)
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 500, verbose = FALSE)
}
info(logger, message = paste("Finished Normalizing Data and Finding Variable Features on Seurat Objects"))
info(logger, message = "==========================================================")

reference_list <- seurat_objects
info(logger, message = "==========================================================")
info(logger, message = paste("Finding Integration Anchors and Integrating Data on Seurat Objects"))
data_anchors <- FindIntegrationAnchors(object.list = reference_list, dims = 1:30)
data_integrated <- IntegrateData(anchorset = data_anchors)
info(logger, message = paste("Finished Finding Integration Anchors and Integrating Data on Seurat Objects"))
info(logger, message = "==========================================================")

DefaultAssay(data_integrated) <- "integrated"

info(logger, message = "==========================================================")
info(logger, message = paste("Scaling Data"))
data_integrated <- ScaleData(data_integrated, verbose = FALSE)
info(logger, message = "Finished Scaling Data")
info(logger, message = "==========================================================")

info(logger, message = "==========================================================")
info(logger, message = paste("Running PCA on Seurat Objects"))
data_integrated <- RunPCA(data_integrated, npcs = args$ncomps, verbose = FALSE)
info(logger, message = paste("Finished running PCA on Seurat Objects"))
info(logger, message = "==========================================================")

info(logger, message = "==========================================================")
info(logger, message = "Running UMAP on Seurat Objects")
data_integrated <- RunUMAP(data_integrated, reduction = "pca", dims = 1:args$ncomps)
info(logger, message = "Finished running UMAP on Seurat Objects")
info(logger, message = "==========================================================")

p1 <- DimPlot(data_integrated, reduction = "umap") + ggplot2::theme(legend.position = "bottom")
p2 <- DimPlot(data_integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf(file = paste(args$output_file_name, ".pdf", sep = ""))
p1 + p2

info(logger, message = "==========================================================")
info(logger, message = "Creating Violin Plot(s) for Seurat Object Features")

if(args$nCount_RNA == TRUE && args$nFeature_RNA == TRUE){
  VlnPlot(data_integrated, c("nCount_RNA", "nFeature_RNA"))
} else if (args$nCount_RNA == TRUE && args$nFeature_RNA == FALSE){
  VlnPlot(data_integrated, c("nCount_RNA"))
} else if (args$nCount_RNA == FALSE && args$nFeature_RNA == TRUE){
  VlnPlot(data_integrated, c("nFeature_RNA"))
}

suppressMessages(suppressWarnings(dev.off()))

info(logger, message = "Created Violin Plot(s) for Seurat Object Features")
info(logger, message = "==========================================================")

info(logger, message = "==========================================================")
info(logger, message = "Saving Integrated Seurat Data as RDS")
saveRDS(data_integrated, file = paste(args$output_file_name, ".rds", sep = ""))
info(logger, message = "Done!")
info(logger, message = "==========================================================")
