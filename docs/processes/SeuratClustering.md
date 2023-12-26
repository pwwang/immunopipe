# SeuratClustering

Cluster all T cells or selected T cells selected by `TCellSelection`.

If `[TCellSelection]` is not set in the configuration, meaning
all cells are T cells, this process will be run on all T cells. Otherwise,
this process will be run on the selected T cells by
[`TCellSelection`](./TCellSelection.md).<br />

See also: [SeuratClusteringOfAllCells](./SeuratClusteringOfAllCells.md).<br />

## Environment Variables

- `ncores` *(`type=int;order=-100`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
    Used in `future::plan(strategy = "multicore", workers = <ncores>)`
    to parallelize some Seurat procedures.<br />
    See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `use_sct` *(`flag;order=-99`)*: *Default: `False`*. <br />
    Whether use SCTransform routine or not
    If `True`, following procedures will be performed in the order:<br />
    * [`SplitObject`](https://satijalab.org/seurat/reference/splitobject).<br />
    * [`SCTransform*`](https://satijalab.org/seurat/reference/sctransform).<br />
    * [`SelectIntegrationFeatures`](https://satijalab.org/seurat/reference/selectintegrationfeatures).<br />
    * [`PrepSCTIntegration`](https://satijalab.org/seurat/reference/prepsctintegration).<br />
    * [`RunPCA*`](https://satijalab.org/seurat/reference/runpca).<br />
    * [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors).<br />
    * [`IntegrateData`](https://satijalab.org/seurat/reference/integratedata).<br />
    * [`RunPCA`](https://satijalab.org/seurat/reference/runpca).<br />
    * [`RunUMAP`](https://satijalab.org/seurat/reference/runumap).<br />
    * [`FindNeighbors`](https://satijalab.org/seurat/reference/findneighbors).<br />
    * [`FindClusters`](https://satijalab.org/seurat/reference/findclusters).<br />
    * `*`: On each sample
    See <https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1>.<br />
    If `False`, fast integration will be performed, using reciprocal PCA (RPCA) and
    following procedures will be performed in the order:<br />
    * [`SplitObject`](https://satijalab.org/seurat/reference/splitobject).<br />
    * [`NormalizeData*`](https://satijalab.org/seurat/reference/normalizedata).<br />
    * [`FindVariableFeatures*`](https://satijalab.org/seurat/reference/findvariablefeatures).<br />
    * [`SelectIntegrationFeatures`](https://satijalab.org/seurat/reference/selectintegrationfeatures).<br />
    * [`ScaleData*`](https://satijalab.org/seurat/reference/scaledata).<br />
    * [`RunPCA*`](https://satijalab.org/seurat/reference/runpca).<br />
    * [`FindIntegrationAnchors`](https://satijalab.org/seurat/reference/findintegrationanchors).<br />
    * [`IntegrateData`](https://satijalab.org/seurat/reference/integratedata).<br />
    * [`ScaleData`](https://satijalab.org/seurat/reference/scaledata).<br />
    * [`RunPCA`](https://satijalab.org/seurat/reference/runpca).<br />
    * [`RunUMAP`](https://satijalab.org/seurat/reference/runumap).<br />
    * [`FindNeighbors`](https://satijalab.org/seurat/reference/findneighbors).<br />
    * [`FindClusters`](https://satijalab.org/seurat/reference/findclusters).<br />
    * `*`: On each sample.<br />
    See <https://satijalab.org/seurat/articles/integration_rpca.html>.<br />
- `SCTransform` *(`ns`)*:
    Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/sctransform>.<br />
- `SelectIntegrationFeatures` *(`ns`)*:
    Arguments for [`SelectIntegrationFeatures()`](https://satijalab.org/seurat/reference/selectintegrationfeatures).<br />
    `object.list` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `nfeatures` *(`type=int`)*: *Default: `3000`*. <br />
        The number of features to select
    - `<more>`:
        See <https://satijalab.org/seurat/reference/selectintegrationfeatures>
- `PrepSCTIntegration` *(`ns`)*:
    Arguments for [`PrepSCTIntegration()`](https://satijalab.org/seurat/reference/prepsctintegration).<br />
    `object.list` and `anchor.features` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/prepsctintegration>
- `NormalizeData` *(`ns`)*:
    Arguments for [`NormalizeData()`](https://satijalab.org/seurat/reference/normalizedata).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/normalizedata>
- `FindVariableFeatures` *(`ns`)*:
    Arguments for [`FindVariableFeatures()`](https://satijalab.org/seurat/reference/findvariablefeatures).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findvariablefeatures>
- `FindIntegrationAnchors` *(`ns`)*:
    Arguments for [`FindIntegrationAnchors()`](https://satijalab.org/seurat/reference/findintegrationanchors).<br />
    `object.list` and `anchor.features` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns for each sample.<br />
    Sample names can also be specified in `reference` instead of indices only.<br />
    `reduction` defaults to `rpca`.<br />
    `normalization.method` defaults to `SCT` if `use_sct` is `True`.<br />
    **If you want to use reference-based integration, you can also set `reference` to a list of sample names, instead of a list of indices.**
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findintegrationanchors>
- `IntegrateData` *(`ns`)*:
    Arguments for [`IntegrateData()`](https://satijalab.org/seurat/reference/integratedata).<br />
    `anchorset` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns for each sample.<br />
    `normalization.method` defaults to `SCT` if `use_sct` is `True`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/integratedata>
- `ScaleData` *(`ns`)*:
    Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata).<br />
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `verbose` *(`flag`)*: *Default: `False`*. <br />
        Whether to print the progress
    - `<more>`:
        See <https://satijalab.org/seurat/reference/scaledata>
- `ScaleData1` *(`ns`)*:
    Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata) that runs on each sample.<br />
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `verbose` *(`flag`)*: *Default: `False`*. <br />
        Whether to print the progress
    - `<more>`:
        See <https://satijalab.org/seurat/reference/scaledata>
- `RunPCA` *(`ns`)*:
    Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca).<br />
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `npcs` *(`type=int`)*:
        The number of PCs to compute.<br />
        For each sample, `npcs` will be no larger than the number of columns - 1.<br />
    - `verbose` *(`flag`)*: *Default: `False`*. <br />
        Whether to print the progress
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runpca>
- `RunPCA1` *(`ns`)*:
    Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca) on each sample.<br />
    `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `npcs` *(`type=int`)*:
        The number of PCs to compute.<br />
        For each sample, `npcs` will be no larger than the number of columns - 1.<br />
    - `verbose` *(`flag`)*: *Default: `False`*. <br />
        Whether to print the progress
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runpca>
- `RunUMAP` *(`ns`)*:
    Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.<br />
    - `dims` *(`type=int`)*: *Default: `30`*. <br />
        The number of PCs to use
    - `reduction`: *Default: `pca`*. <br />
        The reduction to use for UMAP
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runumap>
- `FindNeighbors` *(`ns`)*:
    Arguments for [`FindNeighbors()`](https://satijalab.org/seurat/reference/findneighbors).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findneighbors>
- `FindClusters` *(`ns`)*:
    Arguments for [`FindClusters()`](https://satijalab.org/seurat/reference/findclusters).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `resolution` *(`type=float`)*: *Default: `0.8`*. <br />
        The resolution of the clustering
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findclusters>
- `cache` *(`type=auto`)*: *Default: `True`*. <br />
    Whether to cache the seurat object with cluster information.<br />
    If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.<br />
    The cached seurat object will be saved as `<signature>.cached.RDS` file, where `<signature>` is the signature determined by
    the input and envs of the process.<br />
    See <https://github.com/satijalab/seurat/issues/7849>, <https://github.com/satijalab/seurat/issues/5358> and
    <https://github.com/satijalab/seurat/issues/6748> for more details.<br />
    To not use the cached seurat object, you can either set `cache` to `False` or delete the cached file at
    `.pipen/<Pipeline>/SeuratClustering/0/output/<signature>.cached.RDS`.<br />
    You can also specify the directory to save the cached seurat object by setting `cache` to the directory path.<br />


## Metadata

The metadata of the `Seurat` object will be updated with the cluster
assignments:<br />

![SeuratClustering-metadata](../processes/images/SeuratClustering-metadata.png)
