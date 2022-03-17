# `Seurat.BatchCorrection` Tutorial 📝
---

This tutorial is based off the [**Integration and Label Transfer**](https://satijalab.org/seurat/v3.2/integration.html) Vignette published and developed by the [Satija Lab](https://satijalab.org).

1. `git clone` our `Seurat.BatchCorrection` repository onto your local machine
2. Open up your computer's Terminal and `cd` into the same directory as your `git clone` folder's location. You can also use a graphical interface to navigate to the correct working directory.
4. Once in the correct directory, `cd` into the `/src` folder and run the following command: `Rscript batch_correction.R --input_files demo.txt --use_filenames_for_plots FALSE --ncomps 50 --nCount_RNA TRUE --nFeature_RNA TRUE --output_file_name hnscc_and_mel_bc`
	- `demo.txt` contains the file paths for 2 test datasets within the `/data` directory -- `HNSCC_noribo_small.txt` and `MEL_noribo_small.txt`. These are 500x500 single-cell cancer datasets that will be batch-corrected
	- We utilize several command line flags to provide our users with the freedom to customize their batch correction output
	- To read further about each Command Line Flag / Parameter, take a look at the **technical documentation** for our module here -> [Link](https://github.com/genepattern/Seurat.BatchCorrection/blob/develop/doc.md)
5. You will now have 3 output files: `batch_correction_log.txt`, `hnscc_and_mel_bc.rds`, and `hnscc_and_mel_bc.pdf`
	- The `.txt` file contains a log of each process carried out during the batch correction script's execution
	- The `.rds` file can be used on another one of **GenePattern's** Seurat suite modules, such as the `Seurat.Clustering` module
	- The `.pdf` file contains UMAP plots and Violin plots of the integrated HNSCC and MEL datasets