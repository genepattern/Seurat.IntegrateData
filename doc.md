# Seurat.IntegrateData (v4)
---
**Description**: GenePattern module which implements the batch correction and data integration algorithm derived from the Seurat software package (Seurat version 3.2.0).

**Authors**: Jonathan Zamora and Edwin F. Juárez, UCSD

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)

**Algorithm Version**: Seurat 3.2.0

---

## Summary
---

The `Seurat.IntegrateData` Module integrates (corrects for batch effects) multiple single-cell datasets and identifies shared cell states that are present across different datasets, regardless of their origin. Once the Module integrates these datasets, the returned object will contain a new *Assay* that holds an integrated/batch-corrected expression matrix for all cells. The resultant batch-corrected expression matrix can then be used for downstream analyses and visualizations.


## References
---
"Seurat V3 Paper:" Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM 3rd, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21. doi: [10.1016/j.cell.2019.05.031](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/). Epub 2019 Jun 6. PMID: 31178118; PMCID: PMC6687398.

[Satija Lab](https://satijalab.org)

[Seurat](https://satijalab.org/seurat/)

### Technical Notes
This module runs in the Docker container tagged [`genepattern/seurat-batchcorrection:3.1`](https://hub.docker.com/layers/142293536/genepattern/seurat-batchcorrection/3.1/images/sha256-4b85c00f3d14456a59574acd36b9cfa639f6ce33518abcb9994d12c81a724c09?context=explore).

GitHub repository for this module's code & version: https://github.com/genepattern/Seurat.BatchCorrection/tree/v3

## Parameters
---

| Name | Description |
-------|--------------
| input_files         | Gene expression matrices (columns are cell IDs and rows are genes) stored in a `.txt` or one Seurat Object (`.rds`) file per batch. If your files are in 10x Genomics (Cell Ranger) format or stored in a `.tar`, or `.gz` file, please use the [Seurat.QC](https://cloud.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:genepattern.org:module.analysis:00416) module to pre-process your data and then use the `.txt` or `.rds` output from Seurat.QC in this module. For sample data see [test data file 1](https://datasets.genepattern.org/data/module_support_files/Conos/small_500x500_HNSCC_noribo.txt) and [test data file 2](https://datasets.genepattern.org/data/module_support_files/Conos/small_500x500_MEL_noribo.txt).|
| use_batch_names | (default = TRUE) Map each input file to Batch Numbers beginning with `Batch 1` and going up to `Batch n` where `n` represents the number of input files. When set to FALSE, the batch names will be set to the file names. 
| ncomps            | (default = 50) The number of principal components to be used in the Principal Component Analysis (PCA) for batch correction. |
| nCount_RNA        | (default = TRUE) Whether or not the batch correction script will produce a violin plot of the number of molecules detected within a cell for our single-cell datasets.|
| nFeature_RNA      | (default = TRUE) Whether or not the batch correction script will produce a violin plot of the number of genes detected within each cell of our single-cell datasets.            |
| output_file_name  | (default = 'integrate_data_results') Base name for the `.pdf` and `.rds` output files. |


## Output Files
---

1. `integrate_data_log.txt`
    - This `.txt` file contains a log of each process carried out during the batch correction script's execution.
2. `output_file_name.rds`
    - This `.rds` file can be used on another one of GenePattern's Seurat suite modules, such as the `Seurat.Clustering` module.
    - It contains an integrated / batch-corrected expression matrix for all cells present in the input files.
    - *Note*: This file inherits its name from the `output_file_name` parameter listed above.
3. `output_file_name.pdf`
    - This `.pdf` file contains both a UMAP plot and Violin plot of the integrated Seurat Objects that are derived from gene expression datasets.
    - Additionally, the first page will display a batch mapping table. In short, this table maps input files to numbered batches, starting at `Batch 1` and going up to `Batch n` where n represents the number of input files given to the script.
    - *Note*: This file inherits its name from the `output_file_name` parameter listed above.


## License
---

`Seurat.IntegrateData` is distributed under a modified BSD license available at https://github.com/genepattern/Seurat.IntegrateData/blob/v3/LICENSE.


## Platform Dependencies
---

| Task Type | CPU Type | Operating System | Language |
------------|----------|------------------|----------|
|           |  Any     | Any              | R 4.0.2  |


## Version Comments
---

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 4       | Nov. 29th, 2021 | Updated the naming convention of the module |
| 3       | Mar. 12th, 2021 | Input types now include `.rds`, `tsv`, and `.txt` files. Output shows data before the batch correction (in addition to post batch correction). |
| 2       | Mar. 3rd, 2021 | Introduced updates to parameters, program structure, and output format.|
| 1       | Nov. 4th, 2020 | Initial Release of `Seurat.BatchCorrection`. |
