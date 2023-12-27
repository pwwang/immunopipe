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
- `ScaleData` *(`ns`)*:
    Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata).<br />
    If you want to re-scale the data by regressing to some variables, `Seurat::ScaleData`
    will be called. If nothing is specified, `Seurat::ScaleData` will not be called.<br />
    - `vars-to-regress`:
        The variables to regress on.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/scaledata>
- `SCTransform` *(`ns`)*:
    Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform).<br />
    If you want to re-scale the data by regressing to some variables, `Seurat::SCTransform`
    will be called. If nothing is specified, `Seurat::SCTransform` will not be called.<br />
    - `vars-to-regress`:
        The variables to regress on.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/sctransform>
- `RunUMAP` *(`ns`)*:
    Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.<br />
    - `dims` *(`type=int`)*: *Default: `30`*. <br />
        The number of PCs to use
    - `reduction`:
        The reduction to use for UMAP.<br />
        If not provided, `sobj@misc$integrated_new_reduction` will be used.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runumap>
- `FindNeighbors` *(`ns`)*:
    Arguments for [`FindNeighbors()`](https://satijalab.org/seurat/reference/findneighbors).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    - `reduction`:
        The reduction to use.<br />
        If not provided, `sobj@misc$integrated_new_reduction` will be used.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findneighbors>
- `FindClusters` *(`ns`)*:
    Arguments for [`FindClusters()`](https://satijalab.org/seurat/reference/findclusters).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    The cluster labels will be saved in `seurat_clusters` and prefixed with "c".<br />
    The first cluster will be "c1", instead of "c0".<br />
    - `resolution`: *Default: `0.8`*. <br />
        The resolution of the clustering. You can have multiple resolutions separated by comma.<br />
        The results will be saved in `seurat_clusters_<resolution>`.<br />
        The final resolution will be used to define the clusters at `seurat_clusters`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findclusters>
- `cache` *(`type=auto`)*: *Default: `True`*. <br />
    Whether to cache the seurat object with cluster information.<br />
    If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.<br />
    The cached seurat object will be saved as `<signature>.cached.RDS` file, where `<signature>` is the signature determined by
    the input and envs of the process.<br />
    See -
    * <https://github.com/satijalab/seurat/issues/7849>
    * <https://github.com/satijalab/seurat/issues/5358> and
    * <https://github.com/satijalab/seurat/issues/6748> for more details.<br />
    To not use the cached seurat object, you can either set `cache` to `False` or delete the cached file at
    `<signature>.cached.RDS` in the cache directory.<br />
    If `True`, the cache directory is `.pipen/<Pipeline>/SeuratClustering/0/output/`
    You can also specify customized directory to save the cached seurat object by setting `cache` to the directory path.<br />

## Metadata

The metadata of the `Seurat` object will be updated with the cluster
assignments:<br />

![SeuratClustering-metadata](../processes/images/SeuratClustering-metadata.png)

