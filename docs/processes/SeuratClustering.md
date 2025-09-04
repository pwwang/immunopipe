# SeuratClustering

Cluster all T cells or selected T cells selected by `TCellSelection`.

If `[TCellSelection]` is not set in the configuration, meaning
all cells are T cells, this process will be run on all T cells. Otherwise,
this process will be run on the selected T cells by
[`TCellSelection`](./TCellSelection.md).<br />

See also: [SeuratClusteringOfAllCells](./SeuratClusteringOfAllCells.md).<br />

## Input

- `srtobj`:
    The seurat object loaded by SeuratPreparing

## Output

- `outfile`: *Default: `{{in.srtobj | stem}}.qs`*. <br />
    The seurat object with cluster information at `seurat_clusters`.<br />

## Environment Variables

- `ncores` *(`type=int;order=-100`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
    Used in `future::plan(strategy = "multicore", workers = <ncores>)`
    to parallelize some Seurat procedures.<br />
    See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `RunUMAP` *(`ns`)*:
    Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.<br />
    - `dims` *(`type=int`)*:
        The number of PCs to use
    - `reduction`:
        The reduction to use for UMAP.<br />
        If not provided, `sobj@misc$integrated_new_reduction` will be used.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runumap>
- `RunPCA` *(`ns`)*:
    Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca).<br />
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
    - `resolution` *(`type=auto`)*: *Default: `0.8`*. <br />
        The resolution of the clustering. You can have multiple resolutions as a list or as a string separated by comma.<br />
        Ranges are also supported, for example: `0.1:0.5:0.1` will generate `0.1, 0.2, 0.3, 0.4, 0.5`. The step can be omitted, defaulting to 0.1.<br />
        The results will be saved in `seurat_clusters_<resolution>`.<br />
        The final resolution will be used to define the clusters at `seurat_clusters`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findclusters>
- `cache` *(`type=auto`)*: *Default: `/tmp/m161047`*. <br />
    Where to cache the information at different steps.<br />
    If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.<br />
    Set to `False` to not cache the results.<br />

## Metadata

The metadata of the `Seurat` object will be updated with the cluster
assignments:<br />

![SeuratClustering-metadata](../..//processes/images/SeuratClustering-metadata.png)

