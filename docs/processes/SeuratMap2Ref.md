# SeuratMap2Ref

Map the seurat object to reference

See: <https://satijalab.org/seurat/articles/integration_mapping.html>
and <https://satijalab.org/seurat/articles/multimodal_reference_mapping.html>

## Environment Variables

- `ncores` *(`type=int;order=-100`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
    Used in `future::plan(strategy = "multicore", workers = <ncores>)`
    to parallelize some Seurat procedures.<br />
    See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `use`:
    A column name of metadata from the reference
    (e.g. `celltype.l1`, `celltype.l2`) to transfer to the query as the
    cell types (ident) for downstream analysis. This field is required.<br />
    If you want to transfer multiple columns, you can use
    `envs.MapQuery.refdata`.<br />
- `ident`: *Default: `seurat_clusters`*. <br />
    The name of the ident for query transferred from `envs.use` of the reference.<br />
- `ref`:
    The reference seurat object file.<br />
    Either an RDS file or a h5seurat file that can be loaded by
    `Seurat::LoadH5Seurat()`.<br />
    The file type is determined by the extension. `.rds` or `.RDS` for
    RDS file, `.h5seurat` or `.h5` for h5seurat file.<br />
- `SCTransform` *(`ns`)*:
    Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform)
    - `do-correct-umi` *(`flag`)*: *Default: `False`*. <br />
        Place corrected UMI matrix in assay counts layer?<br />
    - `do-scale` *(`flag`)*: *Default: `False`*. <br />
        Whether to scale residuals to have unit variance?<br />
    - `do-center` *(`flag`)*: *Default: `True`*. <br />
        Whether to center residuals to have mean zero?<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/sctransform>.<br />
        Note that the hyphen (`-`) will be transformed into `.` for the keys.<br />
- `FindTransferAnchors` *(`ns`)*:
    Arguments for [`FindTransferAnchors()`](https://satijalab.org/seurat/reference/findtransferanchors)
    - `normalization-method` *(`choice`)*: *Default: `SCT`*. <br />
        Name of normalization method used.<br />
        - `LogNormalize`:
            Log-normalize the data matrix
        - `SCT`:
            Scale data using the SCTransform method
    - `reference-reduction`: *Default: `spca`*. <br />
        Name of dimensional reduction to use from the reference if running the pcaproject workflow.<br />
        Optionally enables reuse of precomputed reference dimensional reduction.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findtransferanchors>.<br />
        Note that the hyphen (`-`) will be transformed into `.` for the keys.<br />
- `MapQuery` *(`ns`)*:
    Arguments for [`MapQuery()`](https://satijalab.org/seurat/reference/mapquery)
    - `reference-reduction`: *Default: `spca`*. <br />
        Name of reduction to use from the reference for neighbor finding
    - `reduction-model`: *Default: `wnn.umap`*. <br />
        `DimReduc` object that contains the umap model.<br />
    - `refdata` *(`type=json`)*: *Default: `{}`*. <br />
        Extra data to transfer from the reference to the query.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/mapquery>.<br />
        Note that the hyphen (`-`) will be transformed into `.` for the keys.<br />
- `MappingScore` *(`ns`)*:
    Arguments for [`MappingScore()`](https://satijalab.org/seurat/reference/mappingscore)
    - `<more>`:
        See <https://satijalab.org/seurat/reference/mappingscore>.<br />
        Note that the hyphen (`-`) will be transformed into `.` for the keys.<br />
    - `ndim`: *Default: `30`*. <br />

## Metadata

The metadata of the `Seurat` object will be updated with the cluster
assignments (column name determined by `envs.name`):<br />

![SeuratMap2Ref-metadata](../processes/images/SeuratClustering-metadata.png)

