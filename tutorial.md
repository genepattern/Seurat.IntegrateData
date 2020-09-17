# `Seurat.BatchCorrection` Tutorial 📝

1. `git clone` our `Seurat.BatchCorrection` onto your local machine
2. `demo.txt` will reference 2 test datasets -- `HNSCC_noribo_small.txt` and `MEL_noribo_small.txt`. These are 500x500 single-cell cancer datasets that will be batch corrected
3. Open up Terminal and `cd` into the same directory as your `git clone` folder's location
4. Once in the correct directory, run the following command `Rscript batch_correction.R --file_list demo.txt --ncomps 30 --nCount_RNA TRUE --nFeature_RNA TRUE --output_file_name hnscc_and_mel_bc`
	- We utilize several command line flags to provide our users with the freedom to customize their batch correction output
	- `--file_list` is a `.txt` file that specifies the current directory of `.txt` files with single-cell data
	- `ncomps` specifies the number of principal components to be used in the Principal Component Analysis (PCA)
	- `nCount_RNA` is a boolean TRUE / FALSE parameter that will affect whether or not the batch correction script will produce a violin plot of the number of molecules detected within a cell for our single-cell datasets
	- `nFeature_RNA` is a boolean TRUE / FALSE parameter that will affect whether or not the batch correction script will produce a violin plot of the number of genes detected within each cell of our single-cell datasets
	- `output_file_name` simply modifies the name of the `.pdf` and `.rds` output files once our Batch Correction algorithm has finished running
5. You will now have 3 output files: `batch_processing_log.txt`, `hnscc_and_mel_bc.rds`, and `hnscc_and_mel_bc.pdf`
	- The `.rds` file can be used on another one of **GenePattern's** Seurat suite modules, such as the `Seurat.Clustering` module
	- The `.pdf` file contains UMAP plots and Violin plots of the integrated HNSCC and MEL datasets