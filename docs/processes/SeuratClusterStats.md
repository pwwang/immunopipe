# SeuratClusterStats

Statistics of the clustering, including the number/fraction of cells in each cluster, the gene expression values and dimension reduction plots. It's also possible to perform stats on TCR clones/clusters for each T-cell cluster.

## Environment variables

- `stats_defaults` (`ns`): The default parameters for `stats`.
    The parameters from the cases can overwrite the default parameters.
    - `frac` (`flag`): Whether to output the fraction of cells instead of number.
    - `pie` (`flag`): Also output a pie chart?
    - `table` (`flag`): Whether to output a table (in tab-delimited format) and in the report.
    - `ident`: The column name in metadata to use as the identity.
    - `group-by`: The column name in metadata to group the cells.
        Does NOT support for pie charts.
    - `split-by`: The column name in metadata to split the cells into
        different plots.
    - `subset`: An expression to subset the cells, will be passed to
        `dplyr::filter()` on metadata.
    - `devpars` (`ns`): The device parameters for the plots.
        - `res` (`type=int`): The resolution of the plots.
        - `height` (`type=int`): The height of the plots.
        - `width` (`type=int`): The width of the plots.
- `stats` (`type=json`): The number/fraction of cells to plot.
    Keys are the names of the plots and values are the dicts inherited from `env.stats_defaults`.
    Here are some examples -
    ```
    >>> {
    >>>     "nCells_All": {},
    >>>     "nCells_Sample": {"kind": "num", "group-by": "Sample"},
    >>>     "fracCells_Sample": {"kind": "frac", "group-by": "Sample"},
    >>> }
    ```
- `features_defaults` (`ns`): The default parameters for `features`.
    - `features`: The features to plot.
        It can be either a string with comma separated features, a list of features, a file path with `file://` prefix with features
        (one per line), or an integer to use the top N features from [`VariantFeatures(srtobj)`](https://satijalab.org/seurat/reference/VariantFeatures).
    - `ident`: The column name in metadata to use as the identity.
    - `subset`: An expression to subset the cells, will be passed to [`tidyrseurat::filter()`](https://stemangiola.github.io/tidyseurat/reference/filter.html).
    - `devpars` (`ns`): The device parameters for the plots. Does not work for `table`.
        - `res` (`type=int`): The resolution of the plots.
        - `height` (`type=int`): The height of the plots.
        - `width` (`type=int`): The width of the plots.
    - `plus`: The extra elements to add to the `ggplot` object. Does not work for `table`.
    - `group-by`: Group cells in different ways (for example, orig.ident). Works for `ridge`, `vln`, and `dot`.
        It also works for `feature` as `shape.by` being passed to [`Seurat::FeaturePlot`](https://satijalab.org/seurat/reference/featureplot).
    - `split-by`: The column name in metadata to split the cells into different plots.
        It works for `vln`, `feature`, and `dot`.
    - `assay`: The assay to use.
    - `slot`: The slot to use.
    - `section`: The section to put the plot in the report.
        If not specified, the case title will be used.
    - `ncol`: The number of columns for the plots.
    - `kind` (`choice`): The kind of the plot or table.
        - `ridge`: Use [`Seurat::RidgePlot`](https://satijalab.org/seurat/reference/RidgePlot).
        - `ridgeplot`: Same as `ridge`.
        - `vln`: Use [`Seurat::VlnPlot`](https://satijalab.org/seurat/reference/VlnPlot).
        - `vlnplot`: Same as `vln`.
        - `violin`: Same as `vln`.
        - `violinplot`: Same as `vln`.
        - `feature`: Use [`Seurat::FeaturePlot`](https://satijalab.org/seurat/reference/FeaturePlot).
        - `featureplot`: Same as `feature`.
        - `dot`: Use [`Seurat::DotPlot`](https://satijalab.org/seurat/reference/DotPlot).
        - `dotplot`: Same as `dot`.
        - `heatmap`: Use [`Seurat::DoHeatmap`](https://satijalab.org/seurat/reference/DoHeatmap).
            You can specify `average=True` to plot on the average of the expressions.
        - `table`: The table for the features, only gene expressions are supported.
            (supported keys: ident, subset, and features).
- `features` (`type=json`): The plots for features, include gene expressions, and columns from metadata.
    Keys are the titles of the cases and values are the dicts inherited from `env.features_defaults`. It can also have other parameters from
    each Seurat function used by `kind`. Note that for argument name with `.`, you should use `-` instead.
- `dimplots_defaults` (`ns`): The default parameters for `dimplots`.
    - `ident`: The column name in metadata to use as the identity.
    - `group-by`: Same as `ident`. How the points are colored.
    - `split-by`: The column name in metadata to split the cells into different plots.
    - `shape-by`: The column name in metadata to use as the shape.
    - `devpars` (`ns`): The device parameters for the plots.
        - `res` (`type=int`): The resolution of the plots.
        - `height` (`type=int`): The height of the plots.
        - `width` (`type=int`): The width of the plots.
    - `reduction` (`choice`): Which dimensionality reduction to use.
        - `dim`: Use [`Seurat::DimPlot`](https://satijalab.org/seurat/reference/DimPlot).
            First searches for `umap`, then `tsne`, then `pca`.
        - `auto`: Same as `dim`
        - `umap`: Use [`Seurat::UMAPPlot`](https://satijalab.org/seurat/reference/DimPlot).
        - `tsne`: Use [`Seurat::TSNEPlot`](https://satijalab.org/seurat/reference/DimPlot).
        - `pca`: Use [`Seurat::PCAPlot`](https://satijalab.org/seurat/reference/DimPlot).
    - `<more>`: See <https://satijalab.org/seurat/reference/dimplot>
- `dimplots` (`type=json`): The dimensional reduction plots.
    Keys are the titles of the plots and values are the dicts inherited from `env.dimplots_defaults`. It can also have other parameters from
    [`Seurat::DimPlot`](https://satijalab.org/seurat/reference/dimplot).


## Examples

### Number of cells in each cluster

```toml
[SeuratClusterStats.envs.stats]
# suppose you have nothing set in `envs.stats_defaults`
# otherwise, the settings will be inherited here
nCells_All = { }
```

![nCells_All](images/SeuratClusterStats_nCells_All.png){: width="80%" }

### Number of cells in each cluster by groups

```toml
[SeuratClusterStats.envs.stats]
nCells_Sample = { group-by = "Sample" }
```

![nCells_Sample](images/SeuratClusterStats_nCells_Sample.png){: width="80%" }

### Violin plots for the gene expressions

```toml
[SeuratClusterStats.envs.features]
features = "CD4,CD8A"
# Remove the dots in the violin plots
vlnplots = { pt-size = 0, kind = "vln" }
# Don't use the default genes
vlnplots_1 = { features = ["FOXP3", "IL2RA"], pt-size = 0, kind = "vln" }
```

![vlnplots](images/SeuratClusterStats_vlnplots.png){: width="80%" }
![vlnplots_1](images/SeuratClusterStats_vlnplots_1.png){: width="80%" }

### Dimension reduction plot with labels

```toml
[SeuratClusterStats.envs.dimplots.Idents]
label = true
label-box = true
repel = true
```

![dimplots](images/SeuratClusterStats_dimplots.png){: width="80%" }
