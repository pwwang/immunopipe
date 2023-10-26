# SeuratPreparing

Load, prepare and apply QC to data, using `Seurat`

This process will -
- Prepare the seurat object
- Apply QC to the data

See also
- <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1)>
- <https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Create_one_merged_object>

This process will read the scRNA-seq data, based on the information provided by
`SampleInfo`, specifically, the paths specified by the `RNAData` column.<br />
Those paths should be either paths to directoies containing `matrix.mtx`,
`barcodes.tsv` and `features.tsv` files that can be loaded by
[`Seurat::Read10X()`](https://satijalab.org/seurat/reference/read10x),
or paths to `h5` files that can be loaded by
[`Seurat::Read10X_h5()`](https://satijalab.org/seurat/reference/read10x_h5).<br />

Each sample will be loaded individually and then merged into one `Seurat` object, and then perform QC.<br />

In order to perform QC, some additional columns are added to the meta data of the `Seurat` object. They are:<br />

- `precent.mt`: The percentage of mitochondrial genes.<br />
- `percent.ribo`: The percentage of ribosomal genes.<br />
- `precent.hb`: The percentage of hemoglobin genes.<br />
- `percent.plat`: The percentage of platelet genes.<br />

See also [Preparing the input](../preparing-input.md#scRNA-seq-data).<br />

## Environment Variables

- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
    Used in `future::plan(strategy = "multicore", workers = <ncores>)`
    to parallelize some Seurat procedures.<br />
- `cell_qc`:
    Filter expression to filter cells, using
    `tidyrseurat::filter()`.<br />
    Available QC keys include `nFeature_RNA`, `nCount_RNA`,
    `percent.mt`, `percent.ribo`, `percent.hb`, and `percent.plat`.<br />

    /// Tip | Example
    Including the columns added above, all available QC keys include
    `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`, `percent.hb`,
    and `percent.plat`. For example:<br />

    ```toml
    [SeuratPreparing.envs]
    cell_qc = "nFeature_RNA > 200 & percent.mt < 5"
    ```
    will keep cells with more than 200 genes and less than 5%% mitochondrial
    genes.<br />
    ///

- `gene_qc` *(`ns`)*:
    Filter genes. Currently only `min_cells` is supported.<br />
    `gene_qc` is applied after `cell_qc`.<br />
    - `min_cells`: *Default: `3`*. <br />
        The minimum number of cells that a gene must be
        expressed in to be kept.<br />

        /// Tip | Example
        ```toml
        [SeuratPreparing.envs]
        gene_qc = { min_cells = 3 }
        ```
        will keep genes that are expressed in at least 3 cells.<br />
        ///

