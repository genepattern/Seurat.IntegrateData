# Seurat.BatchCorrection (v1)
---
**Description**: GenePattern module which implements the batch correction algorithm derived from the Seurat software package.

**Author**: Jonathan Zamora

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)

**Algorithm Version**: 1

---

## Summary
---

The `Seurat.BatchCorrection` Module aims to integrate / "batch-correct" multiple single-cell datasets and identify shared cell states that are present across different datasets, regardless of their origin. Once the Module integrates these datasets, the returned object will contain a new *Assay* that holds an integrated / batch-corrected expression matrix for all cells. The resultant batch-corrected expression matrix can then be used for downstream analyses and visualizations.


## References
---
[Satija Lab](https://satijalab.org)

[Seurat](https://satijalab.org/seurat/)

### Technical Notes


## Parameters
---

| Name | Description |
-------|--------------
| input_files         | `.txt` file(s) are used as input, and each `.txt` file is composed of gene expression data |
| use_filenames_for_plots | Determines whether `.txt` input file names will be used alongside the batch correction script's UMAP and Violin Plots. By default, this value is `FALSE` since we prefer to map each batch of gene expression data to an ordered batch number beginning with `Batch 1` and going up to `Batch n` where `n` represents the number of input files.
| ncomps            | The number of principal components to be used in the Principal Component Analysis (PCA) in our Batch Correction. By default, our script uses the default value of 50 for the # of Principal Components.                       |
| nCount_RNA        | Boolean TRUE / FALSE parameter that will affect whether or not the batch correction script will produce a violin plot of the number of molecules detected within a cell for our single-cell datasets|
| nFeature_RNA      | Boolean TRUE / FALSE parameter that will affect whether or not the batch correction script will produce a violin plot of the number of genes detected within each cell of our single-cell datasets            |
| output_file_name  | Modifies the name of the `.pdf` and `.rds` output files once our Batch Correction script has finished running                      |


## Output Files
---

1. `batch_correction_log.txt`
    - This `.txt` file contains a log of each process carried out during the batch correction script's execution
2. `output_file_name.rds`
    - This `.rds` file can be used on another one of GenePattern's Seurat suite modules, such as the `Seurat.Clustering` module
    - It contains an integrated / batch-corrected expression matrix for all cells present in the input files
    - *Note*: This file inherits its name from the `output_file_name` parameter listed above.
3. `output_file_name.pdf`
    - This `.pdf` file contains both a UMAP plot and Violin plot of the integrated Seurat Objects that are derived from gene expression datasets
    - Additionally, the first page will display a batch mapping table. In short, this table maps input files to numbered batches, starting at `Batch 1` and going up to `Batch n` where n represents the number of input files given to the script.
    - *Note*: This file inherits its name from the `output_file_name` parameter listed above.


## License
---

`Seurat.BatchCorrection` is distributed under a modified BSD license available at https://github.com/genepattern/Seurat.BatchCorrection/blob/develop/LICENSE


## Platform Dependencies
---

| Task Type | CPU Type | Operating System | Language |
------------|----------|------------------|----------|
|           |  Any     | Any              | R 4.0.2  |


## Version Comments
---

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 1       | TBD          | Initial Release of `Seurat.BatchCorrection` |