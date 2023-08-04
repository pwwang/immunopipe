# SeuratPreparing

This process will read the scRNA-seq data, based on the information provided by `SampleInfo`, specifically, the paths specified by the `RNAData` column. Those paths should be either paths to directoies containing `matrix.mtx`, `barcodes.tsv` and `features.tsv` files that can be loaded by [`Seurat::Read10X()`][1], or paths to `h5` files that can be loaded by [`Seurat::Read10X_h5()`][2].

See also [Preparing the input](../preparing-input.md#scRNA-seq-data).

Each sample will be loaded individually and then merged into one `Seurat` object, and then perform QC.

In order to perform QC, some additional columns are added to the meta data of the `Seurat` object. They are:

- `precent.mt`: The percentage of mitochondrial genes.
- `percent.ribo`: The percentage of ribosomal genes.
- `precent.hb`: The percentage of hemoglobin genes.
- `percent.plat`: The percentage of platelet genes.

## Environment variables

- `cell_qc`: Filter expression to filter cells, using [tidyrseurat::filter()][3].

    !!! tip
        Including the columns added above, all available QC keys include `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`, `percent.hb`, and `percent.plat`. For example:

        ```toml
        [SeuratPreparing.envs]
        cell_qc = "nFeature_RNA > 200 & percent.mt < 5"
        ```
        will keep cells with more than 200 genes and less than 5% mitochondrial genes.

- `gene_qc`: Filter genes. Currently, only `min_cells` is supported.

    !!! example
        ```toml
        [SeuratPreparing.envs]
        gene_qc = { min_cells = 3 }
        ```
        will keep genes that are expressed in at least 3 cells.

See also:
- [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1) and
- [https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Create_one_merged_object](https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Create_one_merged_object)

[1]: https://satijalab.org/seurat/reference/read10x
[2]: https://satijalab.org/seurat/reference/read10x_h5
[3]: https://stemangiola.github.io/tidyseurat/reference/dplyr-methods.html
