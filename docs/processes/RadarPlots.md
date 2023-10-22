# RadarPlots

Radar plots for cell proportion in different clusters.

This process generates the radar plots for the clusters of T cells.<br />
It explores the proportion of cells in different groups (e.g. Tumor vs Blood)
in different T-cell clusters.<br />

## Environment Variables

- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    Mutaters to mutate the metadata of the seurat object. Keys are the column names and values are the expressions to mutate the columns. These new columns will be used to define your cases..<br />
    See also
    [`mutating the metadata`](../configurations.md#mutating-the-metadata).<br />
- `by`:
    Which column to use to separate the cells in different groups.<br />
    `NA`s will be ignored. For example, If you have a column named `Source`
    that marks the source of the cells, and you want to separate the cells
    into `Tumor` and `Blood` groups, you can set `by` to `Source`.<br />
    The there will be two curves in the radar plot, one for `Tumor` and
    one for `Blood`.<br />
- `each`:
    A column with values to separate all cells in different cases
    When specified, the case will be expanded to multiple cases for
    each value in the column.<br />
    If specified, `section` will be ignored, and the case name will
    be used as the section name.<br />
- `order` *(`list`)*:
    The order of the values in `by`. You can also limit
    (filter) the values we have in `by`. For example, if column `Source`
    has values `Tumor`, `Blood`, `Spleen`, and you only want to plot
    `Tumor` and `Blood`, you can set `order` to `["Tumor", "Blood"]`.<br />
    This will also have `Tumor` as the first item in the legend and `Blood`
    as the second item.<br />
- `cluster_col`: *Default: `seurat_clusters`*. <br />
    The column name of the cluster information.<br />
- `cluster_order` *(`list`)*: *Default: `[]`*. <br />
    The order of the clusters.<br />
    You may also use it to filter the clusters. If not given,
    all clusters will be used.<br />
    If the cluster names are integers, use them directly for the order,
    even though a prefix `Cluster` is added on the plot.<br />
- `breaks` *(`list;itype=int`)*: *Default: `[]`*. <br />
    breaks of the radar plots, from 0 to 100.<br />
    If not given, the breaks will be calculated automatically.<br />
- `direction` *(`choice`)*: *Default: `intra-cluster`*. <br />
    Direction to calculate the percentages.<br />
    - `inter-cluster`:
        the percentage of the cells in all groups
        in each cluster (percentage adds up to 1 for each cluster).<br />
    - `intra-cluster`:
        the percentage of the cells in all clusters.<br />
        (percentage adds up to 1 for each group).<br />
- `section`:
    If you want to put multiple cases into a same section
    in the report, you can set this option to the name of the section.<br />
    Only used in the report.<br />
- `devpars` *(`ns`)*:
    The parameters for `png()`
    - `res` *(`type=int`)*: *Default: `100`*. <br />
        The resolution of the plot
    - `height` *(`type=int`)*: *Default: `1200`*. <br />
        The height of the plot
    - `width` *(`type=int`)*: *Default: `1200`*. <br />
        The width of the plot
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    The cases for the multiple radar plots.<br />
    Keys are the names of the cases and values are the arguments for
    the plots (`each`, `by`, `order`, `breaks`, `direction`,
    `cluster_col`, `cluster_order` and `devpars`).<br />
    If not cases are given, a default case will be used, with the
    key `DEFAULT`.<br />
    The keys must be valid string as part of the file name.<br />

## Examples

Let's say we have a metadata like this:<br />

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

With configurations:<br />

```toml
[RadarPlots.envs]
by = "Source"
```

Then we will have a radar plots like this:<br />

![Radar plots](https://pwwang.github.io/immunopipe/processes/images/RadarPlots-default.png)

We can use `each` to separate the cells into different cases:<br />

```toml
[RadarPlots.envs]
by = "Source"
each = "Timepoint"
```

Then we will have two radar plots, one for `Pre` and one for `Post`:<br />

![Radar plots](https://pwwang.github.io/immunopipe/processes/images/RadarPlots-each.png)

Using `cluster_order` to change the order of the clusters and show only the first 3 clusters:<br />

```toml
[RadarPlots.envs]
by = "Source"
cluster_order = ["2", "0", "1"]
breaks = [0, 50, 100]  # also change the breaks
```

![Radar plots cluster_order](https://pwwang.github.io/immunopipe/processes/images/RadarPlots-cluster_order.png)

/// Attention
All the plots used in the examples are just for demonstration purpose. The real plots will have different appearance.<br />
///

