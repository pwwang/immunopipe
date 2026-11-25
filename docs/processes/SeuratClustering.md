# SeuratClustering

Cluster all cells or selected T/B cells selected by `TOrBCellSelection`.

If `[TOrBCellSelection]` is not set in the configuration, meaning
all cells are T/B cells, this process will be run on all T/B cells. Otherwise,
this process will be run on the selected T/B cells by
[`TOrBCellSelection`](./TOrBCellSelection.md).<br />

/// Note

For new reference-based cell type mapping, use
[`CellTypeAnnotation`](./CellTypeAnnotation.md) with
`envs.tool = "mapquery"`.<br />
Existing [`SeuratMap2Ref`](./SeuratMap2Ref.md) configurations remain supported;
see its deprecation notice. If another annotation process runs in the same
pipeline, set `envs.ident` to a distinct name for this process's cluster
column so the results do not overwrite each other.<br />

///

## Input

- `srtobj`:
    The seurat object loaded by SeuratPreparing

## Output

- `outfile`: *Default: `{{in.srtobj | stem}}.qs`*. <br />
    The seurat object with cluster information at `seurat_clusters` or
    the name specified by `envs.ident`

## Environment Variables

- `ncores` *(`type=int;order=-100`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
    Used in `future::plan(strategy = "multicore", workers = <ncores>)`
    to parallelize some Seurat procedures.<br />
    See also: <https://satijalab.org/seurat/articles/future_vignette.html>
- `ident`: *Default: `seurat_clusters`*. <br />
    The name in the metadata to save the cluster labels.<br />
    A shortcut for `envs["FindClusters"]["cluster.name"]`.<br />
- `python`: *Default: `python`*. <br />
    The path to the python executable where `umap-learn` is installed.<br />
    Needed for `RunUMAP` when `umap.method` is set to `umap-learn`.<br />
- `RunUMAP` *(`ns`)*:
    Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.<br />
    You can also specify `features` instead of `dims` to use specific features for UMAP. It can be
    a list with the following fields: `order` (the order of the markers to use for UMAP, e.g. "desc(abs(avg_log2FC))", and
    `n` (the number of total features to use for UMAP, e.g. 30).<br />
    If `features` is a list, it will run `biopipen.utils::RunSeuratDEAnalysis` to get the markers
    for each group, and then select the top `n`/`ngroups` features for each group based on the `order` field.<br />
    If `features` is a numeric value, it will be treated as the `n` field
    in the list above, with the default `order` being "desc(abs(avg_log2FC))".<br />
    - `dims` *(`type=int`)*:
        The number of PCs to use
    - `reduction`:
        The reduction to use for UMAP.<br />
        If not provided, `sobj@misc$integrated_new_reduction` will be used.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runumap>
- `RunTSNE` *(`ns`)*:
    Arguments for [`RunTSNE()`](https://satijalab.org/seurat/reference/runtsne).<br />
    t-SNE is skipped unless this table is non-empty.<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`.<br />
    `features` may be a character vector, a numeric `n`, or a list with `order` and `n`
    (markers selected through `biopipen.utils::RunSeuratDEAnalysis`, as for `RunUMAP`).<br />
    Note that `Seurat::RunTSNE()` has no `assay` argument, so the `features` path always
    reads the default assay.<br />
    `perplexity` defaults to `max(1, min(30, floor((ncells - 1) / 3)))`.<br />
    - `dims` *(`type=int`)*:
        The number of dimensions to use
    - `reduction`:
        The reduction to use. Defaults to `sobj@misc$integrated_new_reduction`, else `"pca"`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runtsne>
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
    The cluster labels will be saved in cluster names and prefixed with "c".<br />
    The first cluster will be "c1", instead of "c0".<br />
    - `resolution` *(`type=auto`)*: *Default: `0.8`*. <br />
        The resolution of the clustering. You can have multiple resolutions as a list or as a string separated by comma.<br />
        Ranges are also supported, for example: `0.1:0.5:0.1` will generate `0.1, 0.2, 0.3, 0.4, 0.5`. The step can be omitted, defaulting to 0.1.<br />
        The results will be saved in `<ident>_<resolution>`.<br />
        The final resolution will be used to define the clusters at `<ident>`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findclusters>
- `cache` *(`type=auto`)*: *Default: `/tmp`*. <br />
    Where to cache the information at different steps.<br />
    If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.<br />
    Set to `False` to not cache the results.<br />

## SeeAlso

- [SeuratClusteringOfAllCells](./SeuratClusteringOfAllCells.md)

## Metadata

The metadata of the `Seurat` object will be updated with the cluster
assignments:<br />

![SeuratClustering-metadata](images/SeuratClustering-metadata.png)

## Description

Clusters the cells to be analysed - all cells, or the T/B cells selected
by `TOrBCellSelection` - using `Seurat` (R), with `FindNeighbors()`,
`FindClusters()` and `RunUMAP()`.<br />

## Base class

`biopipen.ns.scrna.SeuratClustering`

## Deviations

The parameters are passed through to `Seurat` unchanged: this process
overrides no inherited parameter and adds none.<br />

