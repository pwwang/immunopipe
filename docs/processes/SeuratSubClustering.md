# SeuratSubClustering

Sub-cluster the selected T cells.

Find clusters of a subset of cells.<br />

It's unlike [`Seurat::FindSubCluster`], which only finds subclusters of a single
cluster. Instead, it will perform the whole clustering procedure on the subset of
cells. One can use metadata to specify the subset of cells to perform clustering on.<br />

For the subset of cells, the reductions will be re-performed on the subset of cells,
and then the clustering will be performed on the subset of cells. The reduction
will be saved in `sobj@reduction$sub_umap_<casename>` of the original object and the
clustering will be saved in the metadata of the original object using the casename     as the column name.<br />

## Environment Variables

- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata to subset the cells.<br />
    The mutaters will be applied in the order specified.<br />
- `subset`:
    An expression to subset the cells, will be passed to
    [`tidyseurat::filter()`](https://stemangiola.github.io/tidyseurat/reference/filter.html).<br />

- `RunUMAP` *(`ns`)*:
    Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).<br />
    `object` is specified internally as the subset object, and `-` in the key will be replaced with `.`.<br />
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
    - `resolution`: *Default: `0.8`*. <br />
        The resolution of the clustering. You can have multiple resolutions separated by comma.<br />
        The results will be saved in `<casename>_<resolution>`.<br />
        The final resolution will be used to define the clusters at `<casename>`.<br />
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
- `cases` *(`type=json`)*: *Default: `{'subcluster': Diot({})}`*. <br />
    The cases to perform subclustering.<br />
    Keys are the names of the cases and values are the dicts inherited from `envs` except `mutaters` and `cache`.<br />
    If empty, a case with name `subcluster` will be created with default parameters.<br />
