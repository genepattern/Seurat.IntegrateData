# Batch Correction Script
# Written by Jonathan Zamora, UCSD-MesirovLab
# For use in GenePattern

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
# # PARAMETERS for PCA and Violin Plot
# use_filenames_for_plots = FALSE (Default)
# # ====================================
# # PARAMETERS for PCA
# ncomps = 30 (Default)
# # ====================================
# # PARAMETERS for Violin Plot
# nCountRNA: TRUE (Default)
# nFeatureRNA: TRUE (Default)
# # ====================================
# # PARAMETERS for OUTPUT
# output_file_name: "batch_correction_results" (Default)
# # ====================================

# Parse Input Arguments
parser = OptionParser()
# ====================================
parser <- add_option(parser, c("--file_list"), help = "List of files to load. Separate file names with '\ ' ")
parser <- add_option(parser, c("--use_filenames_for_plots"), type = "logical", default = FALSE, help = "Display File Names on UMAP & Violin Plots (default = FALSE)")
# ====================================
# PARAMETERS for PCA
parser <- add_option(parser, c("--ncomps"), type = "integer", default = 30, help = "How many PCA components to use (?) (default = 30)")
# ====================================
# PARAMETERS for Violin Plot
parser <- add_option(parser, c("--nCount_RNA"), type = "logical", default = TRUE, help = "Display nCount_RNA feature on Violin Plot (default = TRUE)")
parser <- add_option(parser, c("--nFeature_RNA"), type = "logical", default = TRUE, help = "Display nFeature_RNA feature on Violin Plot (default = TRUE)")
# ====================================
# PARAMETERS for Output File
parser <- add_option(parser, c("--output_file_name"), default = "batch_correction_results", help = "Output File Name for Batch Correction Analysis (default = 'batch_correction_results')")
# ====================================
info(logger, message = "==========================================================")
args <- parse_args(parser)
info(logger, message = "Parameters used:")
info(logger, message = paste("help:", args$help))
info(logger, message = paste("file_list:", args$file_list))
info(logger, message = paste("use_filenames_for_plots:", args$use_filenames_for_plots))
info(logger, message = paste("ncomps:", args$ncomps))
info(logger, message = paste("nCount_RNA:", args$nCount_RNA))
info(logger, message = paste("nFeature_RNA:", args$nFeature_RNA))
info(logger, message = paste("output_file_name:", args$output_file_name))
info(logger, message = "==========================================================")
# ====================================


# Read file list from command line
# Read the lines of the file and store them in `lines`
con <- file(args$file_list, open = "r")
lines = readLines(con)

# Initialize lists for data and condensed file names
data_list <- list()
condensed_file_names <- list()


info(logger, message = "==========================================================")
info(logger, message = paste("Reading Files contained within:", args$file_list))

for(i in 1:length(lines)){
  data_list[[i]] <- read.table(file = lines[[i]], header = TRUE, sep = "\t", row.names = 1)
  condensed_file_names[[i]] <- tail(strsplit(lines[[i]], "/")[[1]], 1)
}
close(con)

info(logger, message = paste("Read Files contained within:", args$file_list))
info(logger, message = "==========================================================")


info(logger, message = "==========================================================")
info(logger, message = paste("Creating Seurat Objects for files within:", args$file_list))

# Initialize Batch Counter
counter = 1

# Initialize Batch File Name List and Batch Names List
# Also initialize list that holds Seurat Objects
batch_names <- list()
seurat_objects <- list()

# Generate Batch Names in loop
for(i in 1:length(condensed_file_names)){
  batch_names[[i]] <- paste("Batch", counter)
  counter = counter + 1
}

for(i in 1:length(data_list)){
  
  # Default Case for use_filenames_for_plots
  if (args$use_filenames_for_plots == FALSE){
    seurat_objects[[i]] <- CreateSeuratObject(counts = data_list[[i]], project = paste("Batch", counter))
    
    if (length(levels(seurat_objects[[i]])) > 1){
      levels(seurat_objects[[i]]@active.ident) <- rep(batch_names[[i]], length(levels(seurat_objects[[i]])))
    }
    
  } else{
    seurat_objects[[i]] <- CreateSeuratObject(counts = data_list[[i]], project = condensed_file_names[[i]])
    if (length(levels(seurat_objects[[i]])) > 1){
      levels(seurat_objects[[i]]@active.ident) <- rep(condensed_file_names[[i]], length(levels(seurat_objects[[i]])))
    }
  }
}

info(logger, message = paste("Created Seurat Objects for files contained within:", args$file_list))
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


if (args$use_filenames_for_plots == FALSE){
  info(logger, message = "==========================================================")
  info(logger, message = "Creating Dictionary for Seurat Batch Objects")
  pdf(file = paste(args$output_file_name, ".pdf", sep = ""), width = 8.5, height = 11)
  
  # Set up data.frame for Batch Names and their corresponding File Names
  df <- data.frame(do.call(rbind, batch_names))
  colnames(df) <- "Batch_Names"
  df$"File_Names" <- condensed_file_names
  
  
  # These two steps are done in order to correct an error with the data stored in the data.frame
  first.step <- lapply(df, unlist)
  second.step <- as.data.frame(first.step, stringsAsFactors = F)
  
  description = "-----NOTE-----\nThe below table displays a mapping of each input file to a corresponding batch number.\nThis mapping is done in order to simplify the display of the plots on Pg. 2 and Pg. 3."
  plot.new() # Needed in order to use text() function to display description above table
  gridExtra::grid.table(second.step) # Display the Batch Mapping Table
  text(x=0.5, y=1, description, font=2) # Display the description above the Batch Mapping Table
  
  info(logger, message = "Finished Creating Dictionary for Seurat Batch Objects")
  info(logger, message = "==========================================================")
  
} else{ # This will simply generate a pdf file when  args$use_filenames_for_plots == TRUE
  pdf(file = paste(args$output_file_name, ".pdf", sep = ""), width = 8.5, height = 11)
}


info(logger, message = "==========================================================")
info(logger, message = "Plotting UMAP Plots:")
p1 <- DimPlot(data_integrated, reduction = "umap") + ggplot2::theme(legend.position = "bottom")
p2 <- DimPlot(data_integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
p1 + p2
info(logger, message = "Finished Plotting UMAP Plots")
info(logger, message = "==========================================================")


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
