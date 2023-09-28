# RadarPlots

This process generates the radar plots for the clusters of T cells. It explores the proportion of cells in different groups (e.g. Tumor vs Blood) in different T-cell clusters.

## Environment variables

- `mutaters` (`type=json`): Mutaters to mutate the metadata of the
    seurat object. Keys are the column names and values are the
    expressions to mutate the columns. These new columns will be
    used to define your cases.

    See also [`mutating the metadata`][1].

- `by`: Which column to use to separate the cells in different groups.
    `NA`s will be ignored. For example, If you have a column named `Source` that marks the source of the cells, and you want to separate the cells into `Tumor` and `Blood` groups, you can set `by` to `Source`. The there will be two curves in the radar plot, one for `Tumor` and one for `Blood`.

- `order` (`list`): The order of the values in `by`. You can also limit
    (filter) the values we have in `by`. For example, if column `Source` has values `Tumor`, `Blood`, `Spleen`, and you only want to plot `Tumor` and `Blood`, you can set `order` to `["Tumor", "Blood"]`.
    This will also have `Tumor` as the first item in the legend and `Blood` as the second item.

- `each`: A column with values to separate all cells in different cases
    When specified, the case will be expanded to multiple cases for
    each value in the column.
    For example, if column `Timepoint` has values `Pre` and `Post`, then the cells marked with `Pre` will be in one case, and the cells marked with `Post` will be in another case.
    If specified, `section` will be ignored, and the case name will
    be used as the section name.

- `cluster_col`: The column name of the cluster information. Default: `seurat_clusters`.
- `cluster_order`: The order of the clusters to show on the plot in clockwise direction.
- `section`: If you want to put multiple cases into a same section
    in the report, you can set this option to the name of the section.
    Only used in the report.
- `breaks` (`list`;`itype=int`): breaks of the radar plots, 3 integers from 0 to 100.
    If not given, the breaks will be calculated automatically. Literally `[0, max/2, max]`, where `max` is the maximum value of the proportion of cells in the cluster.
- `direction` (`choice`): Direction to calculate the percentages.
    - `inter-cluster`: the percentage of the cells in all groups
        in each cluster (percentage adds up to 1 for each cluster).
    - `intra-cluster`: the percentage of the cells in all clusters.
        (percentage adds up to 1 for each group).
- `devpars` (`ns`): The parameters for `png()` in `R` for plotting.
    - `res` (`type=int`): The resolution of the plot, default: `100`.
    - `height` (`type=int`): The height of the plot, default: `1000`.
    - `width` (`type=int`): The width of the plot, default: `1000`.
- `cases` (`type=json`): The cases for the multiple radar plots.
    Keys are the names of the cases and values are the arguments for
    the plots (`each`, `by`, `order`, `breaks`, `direction`,
    `cluster_col`, `cluster_order` and `devpars`).
    If not cases are given, a default case will be used, with the
    key `DEFAULT`.
    The keys must be valid string as part of the file name.

## Examples

Let's say we have a metadata like this:

| Cell | Source | Timepoint | seurat_clusters |
| ---- | ------ | --------- | --------------- |
| A    | Blood  | Pre       | 0               |
| B    | Blood  | Pre       | 0               |
| C    | Blood  | Post      | 1               |
| D    | Blood  | Post      | 1               |
| E    | Tumor  | Pre       | 2               |
| F    | Tumor  | Pre       | 2               |
| G    | Tumor  | Post      | 3               |
| H    | Tumor  | Post      | 3               |

With configurations:

```toml
[RadarPlots.envs]
by = "Source"
```

Then we will have a radar plots like this:

![Radar plots](images/RadarPlots-default.png)

We can use `each` to separate the cells into different cases:

```toml
[RadarPlots.envs]
by = "Source"
each = "Timepoint"
```

Then we will have two radar plots, one for `Pre` and one for `Post`:

![Radar plots](images/RadarPlots-each.png)

Using `cluster_order` to change the order of the clusters and show only the first 3 clusters:

```toml
[RadarPlots.envs]
by = "Source"
cluster_order = ["2", "0", "1"]
breaks = [0, 50, 100]  # also change the breaks
```

![Radar plots cluster_order](images/RadarPlots-cluster_order.png)


/// Attention
All the plots used in the examples are just for demonstration purpose. The real plots will have different appearance.
///

[1]: ../configurations.md#mutating-the-metadata