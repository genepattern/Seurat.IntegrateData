# Integrate Data Script
# Written by Jonathan Zamora, UCSD-MesirovLab
# For use in GenePattern


suppressMessages(suppressWarnings(library("log4r")))
logfile <- "integrate_data_log.txt"
console_appender <- console_appender(layout = default_log_layout())
file_appender <- file_appender(logfile, append = TRUE, layout = default_log_layout())
logger <- log4r::logger(threshold = "INFO", appenders = list(console_appender, file_appender))


info(logger, message = "==========================================================")
info(logger, message = "Loading libraries: Seurat, ggplot2, optparse, log4r")
suppressMessages(suppressWarnings(library("Seurat")))
suppressMessages(suppressWarnings(library("optparse")))
suppressMessages(suppressWarnings(library("ggplot2")))
info(logger, message = "Loaded libraries: Seurat, ggplot2, optparse, log4r")
info(logger, message = "==========================================================")

# # ====================================
# # PARAMETERS for INPUT, PCA, and Violin Plot
# use_batch_names = TRUE (Default)
# # ====================================
# # PARAMETERS for PCA
# ncomps = 50 (Default)
# # ====================================
# # PARAMETERS for Violin Plot
# nCount_RNA: TRUE (Default)
# nFeature_RNA: TRUE (Default)
# # ====================================
# # PARAMETERS for OUTPUT
# output_file_name: "integrate_data_results" (Default)
# # ====================================

# Parameter Descriptions can be found at: https://github.com/genepattern/Seurat.IntegrateData/blob/develop/doc.md
# ====================================
# Parse Input Arguments
parser = OptionParser()
# ====================================
parser <- add_option(parser, c("--input_files"), help="List of files to load. Separate file names with '\ ' ")
parser <- add_option(parser, c("--use_batch_names"), type="logical", default=TRUE, help="Maps each input file to Batch Numbers (default = TRUE)")
# ====================================
# PARAMETERS for PCA
parser <- add_option(parser, c("--ncomps"), type="integer", default=50, help = "How many PCA components to use (default = 50)")
# ====================================
# PARAMETERS for Violin Plot
parser <- add_option(parser, c("--nCount_RNA"), type="logical", default=TRUE, help="Display nCount_RNA feature on Violin Plot (default = TRUE)")
parser <- add_option(parser, c("--nFeature_RNA"), type="logical", default=TRUE, help="Display nFeature_RNA feature on Violin Plot (default = TRUE)")
# ====================================
# PARAMETERS for Output File
parser <- add_option(parser, c("--output_file_name"), default="batch_correction_results", help="Output File Name for Data Integration Analysis (default = 'integrate_data_results')")
# ====================================


info(logger, message = "==========================================================")
info(logger, message = "Parameters used:")
args <- parse_args(parser)
info(logger, message = paste("help:", args$help))
info(logger, message = paste("input_files:", args$input_files))
info(logger, message = paste("use_batch_names:", args$use_batch_names))
info(logger, message = paste("ncomps:", args$ncomps))
info(logger, message = paste("nCount_RNA:", args$nCount_RNA))
info(logger, message = paste("nFeature_RNA:", args$nFeature_RNA))
info(logger, message = paste("output_file_name:", args$output_file_name))
info(logger, message = "==========================================================")


con <- file(args$input_files, open = "r") # Open "read" connection to input files
lines = readLines(con) # Read lines of input files from open "read" connection

# Initialize lists for data, data labels, and batch names
data <- list()
data_labels <- list()
batch_names <- list()
counter = 1 # Initialize Batch Counter

info(logger, message = "==========================================================")
info(logger, message = paste("Reading Files contained within:", args$input_files))

# Run this code by default
for(i in 1:length(lines)){
  
  if (grepl("\\.txt$", lines[[i]]) || grepl("\\.tsv$", lines[[i]]))
  {
    data[[i]] <- read.delim(file = lines[[i]], row.names=1) # creates data frame for each input file
    data_labels[[i]] <- tail(strsplit(lines[[i]], "/")[[1]], 1) # stores corresponding data labels
    
    if (args$use_batch_names == TRUE ) {
      batch_names[[i]] <- paste("Batch", counter) # store batch names in list
    }
    
    else if (args$use_batch_names == FALSE) {
      batch_names[[i]] <- tail(strsplit(lines[[i]], "/")[[1]], 1) # stores data labels for each batch name
    }
    counter = counter + 1 # increment batch counter
  } else if (grepl("\\.rds$", lines[[i]])) {
    data[[i]] <- readRDS(file = lines[[i]]) # creates data frame for each input file
    data_labels[[i]] <- tail(strsplit(lines[[i]], "/")[[1]], 1) # stores corresponding data labels
    
    if (args$use_batch_names == TRUE ) {
      batch_names[[i]] <- paste("Batch", counter) # store batch names in list
    }
    
    else if (args$use_batch_names == FALSE) {
      batch_names[[i]] <- tail(strsplit(lines[[i]], "/")[[1]], 1) # stores data labels for each batch name
    }
    counter = counter + 1 # increment batch counter
  } else {
    info(logger, message = paste("ERROR: The current file", lines[[i]], "is not in a valid format. Please use .txt, .tsv, or .rds file inputs only"))
    info(logger, message = paste("If you have .tar.gz or .zip files, consult Seurat.QC [https://cloud.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:genepattern.org:module.analysis:00416:2] to convert those files to .txt or .rds format"))
    info(logger, message = paste("END OF PROGRAM"))
    quit(save="no")
  }
}

close(con) # Close "read" connection for input files; No longer needed

info(logger, message = paste("Finished Reading Files contained within:", args$input_files))
info(logger, message = "==========================================================")


info(logger, message = "==========================================================")
info(logger, message = paste("Creating Seurat Objects for files within:", args$input_files))

seurat_objects <- list() # Initialize Seurat Object list

for(i in 1:length(data))
{
  if (is(data[[i]], "Seurat") == TRUE && args$use_batch_names == TRUE) {
    seurat_objects[[i]] <- data[[i]]
    Project(seurat_objects[[i]]) <- batch_names[[i]]
    # Associates all levels in Seurat Object with Current Batch
    levels(seurat_objects[[i]]@active.ident) <- rep(batch_names[[i]], length(levels(seurat_objects[[i]])))
  }
  
  else if (is(data[[i]], "Seurat") == TRUE && args$use_batch_names == FALSE) {
    seurat_objects[[i]] <- data[[i]]
    Project(seurat_objects[[i]]) <- data_labels[[i]]
    levels(seurat_objects[[i]]@active.ident) <- rep(data_labels[[i]], length(levels(seurat_objects[[i]])))
  }
  
  else {
    if (args$use_batch_names == TRUE) {
      seurat_objects[[i]] <- CreateSeuratObject(counts = data[[i]], project = batch_names[[i]])
      
      # Associates all levels in Seurat Object with Current Batch
      levels(seurat_objects[[i]]@active.ident) <- rep(batch_names[[i]], length(levels(seurat_objects[[i]])))
      
    } else {
      seurat_objects[[i]] <- CreateSeuratObject(counts = data[[i]], project = data_labels[[i]])
      levels(seurat_objects[[i]]@active.ident) <- rep(data_labels[[i]], length(levels(seurat_objects[[i]])))
    }
  }
}

batch_names <- append(batch_names, paste("Batch", (length(batch_names) + 1))) # Add batch number for merged seurat objects
merged_seurat <- merge(x=seurat_objects[[1]], y=seurat_objects[2:length(seurat_objects)], project = batch_names[[length(batch_names)]]) # merge (but dont correct) seurat objects
seurat_objects <- append(seurat_objects, merged_seurat) # append merged seurat object to end of seurat object list

info(logger, message = paste("Created Seurat Objects for files contained within:", args$input_files))
info(logger, message = "==========================================================")


info(logger, message = "==========================================================")
info(logger, message = paste("Normalizing Data and Finding Variable Features on Seurat Objects"))

for(i in 1:length(seurat_objects)){
  seurat_objects[[i]] <- NormalizeData(seurat_objects[[i]], verbose = FALSE)
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 500, verbose = FALSE)
}

info(logger, message = paste("Finished Normalizing Data and Finding Variable Features on Seurat Objects"))
info(logger, message = "==========================================================")


reference_list <- seurat_objects[1:length(seurat_objects)-1]


info(logger, message = "==========================================================")
info(logger, message = paste("Finding Integration Anchors and Integrating Data on Seurat Objects"))
data_anchors <- FindIntegrationAnchors(object.list = reference_list, dims=1:args$ncomps)
data_integrated <- IntegrateData(anchorset = data_anchors)
info(logger, message = paste("Finished Finding Integration Anchors and Integrating Data on Seurat Objects"))
info(logger, message = "==========================================================")


DefaultAssay(data_integrated) <- "integrated" # Set name of default assay as "integrated"


info(logger, message = "==========================================================")
info(logger, message = paste("Scaling Data"))
data_integrated <- ScaleData(data_integrated, verbose = FALSE) # Scale and Center Features in dataset(s)
info(logger, message = "Finished Scaling Data")
info(logger, message = "==========================================================")


info(logger, message = "==========================================================")
info(logger, message = paste("Running PCA on Seurat Objects"))
data_integrated <- RunPCA(data_integrated, npcs = args$ncomps, verbose = FALSE, seed.use = 42) # Dimensionality Reduction
info(logger, message = paste("Finished running PCA on Seurat Objects"))
info(logger, message = "==========================================================")


info(logger, message = "==========================================================")
info(logger, message = "Running UMAP on Seurat Objects")
data_integrated <- RunUMAP(data_integrated, reduction.use = "pca", dims = 1:args$ncomps, seed.use = 42) # Another Dimensional Reduction Technique
info(logger, message = "Finished running UMAP on Seurat Objects")
info(logger, message = "==========================================================")


info(logger, message = "==========================================================")
info(logger, message = "Creating Dictionary for Seurat Batch Objects")
pdf(file = paste(args$output_file_name, ".pdf", sep = ""), width = 8.5, height = 11)

# Set up data.frame for Batch Names and their corresponding File Names
df <- data.frame(do.call(rbind, batch_names[1:length(batch_names)-1]))
colnames(df) <- "Batch_Names"
df$"File_Names" <- data_labels


# These two steps are executed in order to correct an error with how the data is originally stored in the data.frame
first.step <- lapply(df, unlist)
second.step <- as.data.frame(first.step, stringsAsFactors = F)

description = "NOTE:\nThe below table displays a mapping of each input file to a corresponding batch number.\nThis mapping is completed in order to simplify the display of the plots on Pg. 2 and Pg. 3."
plot.new() # Needed in order to use text() function to display description above table
gridExtra::grid.table(second.step) # Display the Batch Mapping Table
text(x=0, y=1, description, font=2, pos = 4, cex = 0.9) # Display the description above the Batch Mapping Table

info(logger, message = "Finished Creating Dictionary for Seurat Batch Objects")
info(logger, message = "==========================================================")



info(logger, message = "==========================================================")
info(logger, message = "Plotting UMAP Plots:")

seurat_objects[[length(seurat_objects)]] <- ScaleData(seurat_objects[[length(seurat_objects)]], verbose = FALSE)
seurat_objects[[length(seurat_objects)]] <- RunPCA(seurat_objects[[length(seurat_objects)]], npcs = args$ncomps, verbose = FALSE, seed.use = 42)
seurat_objects[[length(seurat_objects)]] <- RunUMAP(seurat_objects[[length(seurat_objects)]], reduction.use = "pca", dims = 1:args$ncomps, seed.use = 42)

p1 <- DimPlot(seurat_objects[[length(seurat_objects)]], reduction = "umap") + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + ggtitle("Non-Batch-Corrected Data")
p1

p2 <- DimPlot(data_integrated, reduction = "umap") + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + ggtitle("Batch-Corrected Data")
p2
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