# Seurat.BatchCorrection (v1)
---
**Description**: GenePattern module which implements the batch correction algorithm within Seurat.

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
| file_list         | `.txt` file input that contains a list of one or more file paths that point to `.txt` files that are composed of gene expression data |
| ncomps            | The number of principal components to be used in the Principal Component Analysis (PCA) in our Batch Correction                       |
| nCount_RNA        | Boolean TRUE / FALSE parameter that will affect whether or not the batch correction script will produce a violin plot of the number of molecules detected within a cell for our single-cell datasets|
| nFeature_RNA      | Boolean TRUE / FALSE parameter that will affect whether or not the batch correction script will produce a violin plot of the number of genes detected within each cell of our single-cell datasets            |
| output_file_name  | Modifies the name of the `.pdf` and `.rds` output files once our Batch Correction algorithm has finished running                      |


## Output Files
---

1. `batch_processing_log.txt`
    - The `.txt` file contains a log of each process carried out during the batch correction script's execution
2. `<your_output_file_name>.rds`
    - The `.rds` file can be used on another one of GenePattern's Seurat suite modules, such as the `Seurat.Clustering` module
3. `<your_output_file_name>.pdf`
    - The `.pdf` file contains UMAP plots and Violin plots of the integrated Seurat Objects that are derived from gene expression datasets


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