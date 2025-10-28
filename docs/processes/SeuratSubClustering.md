# SeuratSubClustering

Sub-clustering for all or selected T/B cells.

Find clusters of a subset of cells.<br />

It's unlike [`Seurat::FindSubCluster`], which only finds subclusters of a single
cluster. Instead, it will perform the whole clustering procedure on the subset of
cells. One can use metadata to specify the subset of cells to perform clustering on.<br />

For the subset of cells, the reductions will be re-performed on the subset of cells,
and then the clustering will be performed on the subset of cells. The reduction
will be saved in `object@reduction$<casename>.<reduction>` of the original object and the
clustering will be saved in the metadata of the original object using the casename
as the column name.<br />

## Input

- `srtobj`:
    The seurat object in RDS or qs/qs2 format.<br />

## Output

- `outfile`: *Default: `{{in.srtobj | stem}}.qs`*. <br />
    The seurat object with the subclustering information in qs/qs2 format.<br />

## Environment Variables

- `ncores` *(`type=int;order=-100`)*: *Default: `1`*. <br />
    Number of cores to use.<br />
    Used in `future::plan(strategy = "multicore", workers = <ncores>)`
    to parallelize some Seurat procedures.<br />
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata to subset the cells.<br />
    The mutaters will be applied in the order specified.<br />
- `subset`:
    An expression to subset the cells, will be passed to
    [`tidyseurat::filter()`](https://stemangiola.github.io/tidyseurat/reference/filter.html).<br />
- `RunPCA` *(`ns`)*:
    Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca).<br />
    `object` is specified internally as the subset object, and `-` in the key will be replaced with `.`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/runpca>
- `RunUMAP` *(`ns`)*:
    Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).<br />
    `object` is specified internally as the subset object, and `-` in the key will be replaced with `.`.<br />
    `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.<br />
    - `dims` *(`type=int`)*:
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
        If not provided, `object@misc$integrated_new_reduction` will be used.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findneighbors>
- `FindClusters` *(`ns`)*:
    Arguments for [`FindClusters()`](https://satijalab.org/seurat/reference/findclusters).<br />
    `object` is specified internally, and `-` in the key will be replaced with `.`.<br />
    The cluster labels will be prefixed with "s". The first cluster will be "s1", instead of "s0".<br />
    - `resolution` *(`type=auto`)*: *Default: `0.8`*. <br />
        The resolution of the clustering. You can have multiple resolutions as a list or as a string separated by comma.<br />
        Ranges are also supported, for example: `0.1:0.5:0.1` will generate `0.1, 0.2, 0.3, 0.4, 0.5`. The step can be omitted, defaulting to 0.1.<br />
        The results will be saved in `<casename>_<resolution>`.<br />
        The final resolution will be used to define the clusters at `<casename>`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/findclusters>
- `cache` *(`type=auto`)*: *Default: `/tmp`*. <br />
    Whether to cache the results.<br />
    If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.<br />
    Set to `False` to not cache the results.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    The cases to perform subclustering.<br />
    Keys are the names of the cases and values are the dicts inherited from `envs` except `mutaters` and `cache`.<br />
    If empty, a case with name `subcluster` will be created with default parameters.<br />
    The case name will be passed to `biopipen.utils::SeuratSubCluster()` as `name`.<br />
    It will be used as the prefix for the reduction name, keys and cluster names.<br />
    For reduction keys, it will be `toupper(<name>)` + "PC_" and `toupper(<name>)` + "UMAP_".<br />
    For cluster names, it will be `<name>` + "." + resolution.<br />
    And the final cluster name will be `<name>`.<br />
    Note that the `name` should be alphanumeric and anything other than alphanumeric will be removed.<br />

## Metadata

The metadata of the `Seurat` object will be updated with the sub-clusters
specified by names (keys) of `envs.cases`:<br />

![SeuratSubClustering-metadata](images/SeuratSubClustering-metadata.png)

