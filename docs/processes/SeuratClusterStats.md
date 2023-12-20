# SeuratClusterStats

Statistics of the clustering.

Including the number/fraction of cells in each cluster, the gene expression values
and dimension reduction plots. It's also possible to perform stats on
TCR clones/clusters or other metadata for each T-cell cluster.<br />

## Environment Variables

- `stats_defaults` *(`ns`)*:
    The default parameters for `stats`.<br />
    The parameters from the cases can overwrite the default parameters.<br />
    - `frac` *(`flag`)*: *Default: `False`*. <br />
        Whether to output the fraction of cells instead of number.<br />
    - `pie` *(`flag`)*: *Default: `False`*. <br />
        Also output a pie chart?<br />
    - `table` *(`flag`)*: *Default: `False`*. <br />
        Whether to output a table (in tab-delimited format) and in the report.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
    - `group-by`:
        The column name in metadata to group the cells.<br />
        Does NOT support for pie charts.<br />
    - `split-by`:
        The column name in metadata to split the cells into
        different plots.<br />
    - `subset`:
        An expression to subset the cells, will be passed to
        `dplyr::filter()` on metadata.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `800`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `1000`*. <br />
            The width of the plots.<br />
- `stats` *(`type=json`)*: *Default: `{'Number of cells in each cluster': Diot({'pie': True}), 'Number of cells in each cluster by Sample': Diot({'group-by': 'Sample', 'table': True, 'frac': True})}`*. <br />
    The number/fraction of cells to plot.<br />
    Keys are the names of the plots and values are the dicts inherited from `env.stats_defaults`.<br />
    Here are some examples -

    ```python
    {
        "nCells_All": {},
        "nCells_Sample": {"kind": "num", "group-by": "Sample"},
        "fracCells_Sample": {"kind": "frac", "group-by": "Sample"},
    }
    ```

- `ngenes_defaults` *(`ns`)*:
    The default parameters for `ngenes`.<br />
    The default parameters to plot the number of genes expressed in each cell.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
    - `group-by`:
        The column name in metadata to group the cells.<br />
        Dodge position will be used to separate the groups.<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `tidyrseurat::filter()`.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `800`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `1000`*. <br />
            The width of the plots.<br />
- `ngenes` *(`type=json`)*: *Default: `{'Number of genes expressed in each cluster': Diot({})}`*. <br />
    The number of genes expressed in each cell.<br />
    Keys are the names of the plots and values are the dicts inherited from `env.ngenes_defaults`.<br />
- `features_defaults` *(`ns`)*:
    The default parameters for `features`.<br />
    - `features`:
        The features to plot.<br />
        It can be either a string with comma separated features, a list of features, a file path with `file://` prefix with features
        (one per line), or an integer to use the top N features from `VariantFeatures(srtobj)`.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `tidyrseurat::filter()`.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots. Does not work for `table`.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*:
            The height of the plots.<br />
        - `width` *(`type=int`)*:
            The width of the plots.<br />
    - `plus`:
        The extra elements to add to the `ggplot` object. Does not work for `table`.<br />
    - `group-by`:
        Group cells in different ways (for example, orig.ident). Works for `ridge`, `vln`, and `dot`.<br />
        It also works for `feature` as `shape.by` being passed to [`Seurat::FeaturePlot`](https://satijalab.org/seurat/reference/featureplot).<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
        It works for `vln`, `feature`, and `dot`.<br />
    - `assay`:
        The assay to use.<br />
    - `slot`:
        The slot to use.<br />
    - `section`:
        The section to put the plot in the report.<br />
        If not specified, the case title will be used.<br />
    - `ncol` *(`type=int`)*: *Default: `2`*. <br />
        The number of columns for the plots.<br />
    - `kind` *(`choice`)*:
        The kind of the plot or table.<br />
        - `ridge`:
            Use `Seurat::RidgePlot`.<br />
        - `ridgeplot`:
            Same as `ridge`.<br />
        - `vln`:
            Use `Seurat::VlnPlot`.<br />
        - `vlnplot`:
            Same as `vln`.<br />
        - `violin`:
            Same as `vln`.<br />
        - `violinplot`:
            Same as `vln`.<br />
        - `feature`:
            Use `Seurat::FeaturePlot`.<br />
        - `featureplot`:
            Same as `feature`.<br />
        - `dot`:
            Use `Seurat::DotPlot`.<br />
        - `dotplot`:
            Same as `dot`.<br />
        - `heatmap`:
            Use `Seurat::DoHeatmap`.<br />
            You can specify `average=True` to plot on the average of the expressions.<br />
        - `table`:
            The table for the features, only gene expressions are supported.<br />
            (supported keys: ident, subset, and features).<br />
- `features` *(`type=json`)*: *Default: `{}`*. <br />
    The plots for features, include gene expressions, and columns from metadata.<br />
    Keys are the titles of the cases and values are the dicts inherited from `env.features_defaults`. It can also have other parameters from
    each Seurat function used by `kind`. Note that for argument name with `.`, you should use `-` instead.<br />
- `dimplots_defaults` *(`ns`)*:
    The default parameters for `dimplots`.<br />
    - `ident`: *Default: `seurat_clusters`*. <br />
        The column name in metadata to use as the identity.<br />
        Ignored if `group-by` is specified.<br />
    - `group-by`:
        Same as `ident`. How the points are colored.<br />
    - `split-by`:
        The column name in metadata to split the cells into different plots.<br />
    - `shape-by`:
        The column name in metadata to use as the shape.<br />
    - `subset`:
        An expression to subset the cells, will be passed to `tidyrseurat::filter()`.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plots.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plots.<br />
        - `height` *(`type=int`)*: *Default: `800`*. <br />
            The height of the plots.<br />
        - `width` *(`type=int`)*: *Default: `1000`*. <br />
            The width of the plots.<br />
    - `reduction` *(`choice`)*: *Default: `dim`*. <br />
        Which dimensionality reduction to use.<br />
        - `dim`:
            Use `Seurat::DimPlot`.<br />
            First searches for `umap`, then `tsne`, then `pca`.<br />
        - `auto`:
            Same as `dim`
        - `umap`:
            Use `Seurat::UMAPPlot`.<br />
        - `tsne`:
            Use `Seurat::TSNEPlot`.<br />
        - `pca`:
            Use `Seurat::PCAPlot`.<br />
    - `<more>`:
        See <https://satijalab.org/seurat/reference/dimplot>
- `dimplots` *(`type=json`)*: *Default: `{'Dimensional reduction plot': Diot({'label': True, 'label-box': True, 'repel': True}), 'TCR presence': Diot({'ident': 'TCR_Presence', 'order': 'TCR_absent', 'cols': ['#FF000066', 'gray']})}`*. <br />
    The dimensional reduction plots.<br />
    Keys are the titles of the plots and values are the dicts inherited from `env.dimplots_defaults`. It can also have other parameters from
    [`Seurat::DimPlot`](https://satijalab.org/seurat/reference/dimplot).<br />

## Examples

### Number of cells in each cluster

```toml
[SeuratClusterStats.envs.stats]
# suppose you have nothing set in `envs.stats_defaults`
# otherwise, the settings will be inherited here
nCells_All = { }
```

![nCells_All](../processes/images/SeuratClusterStats_nCells_All.png){: width="80%" }

### Number of cells in each cluster by groups

```toml
[SeuratClusterStats.envs.stats]
nCells_Sample = { group-by = "Sample" }
```

![nCells_Sample](../processes/images/SeuratClusterStats_nCells_Sample.png){: width="80%" }

### Violin plots for the gene expressions

```toml
[SeuratClusterStats.envs.features]
features = "CD4,CD8A"
# Remove the dots in the violin plots
vlnplots = { pt-size = 0, kind = "vln" }
# Don't use the default genes
vlnplots_1 = { features = ["FOXP3", "IL2RA"], pt-size = 0, kind = "vln" }
```

![vlnplots](../processes/images/SeuratClusterStats_vlnplots.png){: width="80%" }
![vlnplots_1](../processes/images/SeuratClusterStats_vlnplots_1.png){: width="80%" }

### Dimension reduction plot with labels

```toml
[SeuratClusterStats.envs.dimplots.Idents]
label = true
label-box = true
repel = true
```

![dimplots](../processes/images/SeuratClusterStats_dimplots.png){: width="80%" }

