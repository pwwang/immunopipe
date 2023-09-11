# CellsDistribution

This process will generate a chart with a set of pie charts showing the distribution of cells from different cell groups (i.e. TCR clones) in different seurat clusters (cell sub-types).

## Environment variables

- `mutaters` (`type=json`): The mutaters to mutate the metadata
    Keys are the names of the mutaters and values are the R expressions
    passed by `dplyr::mutate()` to mutate the metadata.
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`, that can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).
    For example, you can use `{"Patient1_Tumor_Collapsed_Clones": "expanded(Source, 'Tumor', subset = Patent == 'Patient1')"}`
    to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
    with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1. The values in this columns for other clones will be `NA`.
    The `expanded` and `contracted` functions take 3 arguments:
    * `df`: The metadata data frame. You can use the `.` to refer to it.
    * `group-by`: The column name in metadata to group the cells.
    * `idents`: The first group or both groups of cells to compare (value in `group-by` column). If only the first group is given, the rest of the cells (with non-NA in `group-by` column) will be used as the second group.
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`)
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df %>% mutate(expanded = expanded(...))`.
    * `order`: The order of the returned ids. It could be `sum` or `diff`, which is the sum or diff of the `compare` between idents.
        Two kinds of modifiers can be added, including `desc` and `abs`.
        For example, `sum,desc` means the sum of `compare` between idents in descending order.
        Default is `diff,abs,desc`. It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
        ids will be in the same order as in `df`.
    Note that the numeric column should be the same for all cells in the same group. This will not be checked (only the first value is used).
- `group_by`: The column name in metadata to group the cells for the columns of the plot.
- `group_order` (`list`): The order of the groups (columns) to show on the plot
- `cells_by`: The column name in metadata to group the cells for the rows of the plot.
- `cells_order` (`list`): The order of the cells (rows) to show on the plot
- `cells_orderby`: An expression passed to `dplyr::arrange()` to order the cells (rows) of the plot.
    Only works when `cells-order` is not specified.
    4 extra columns were added to the metadata for ordering the rows in the plot:
    * `CloneSize`: The size (number of cells) of clones (identified by `cells_by`)
    * `CloneGroupSize`: The clone size in each group (identified by `group_by`)
    * `CloneClusterSize`: The clone size in each cluster (identified by `seurat_clusters`)
    * `CloneGroupClusterSize`: The clone size in each group and cluster (identified by `group_by` and `seurat_clusters`)
- `cells_n` (`type=int`): The max number of groups to show for each cell group identity (row).
    Ignored if `cells_order` is specified.
- `devpars` (`ns`): The device parameters for the plots.
    - `res` (`type=int`): The resolution of the plots
    - `height` (`type=int`): The height of the plots
    - `width` (`type=int`): The width of the plots
- `each`: The column name in metadata to separate the cells into different plots.
- `section`: The section to show in the report. This allows different cases to be put in the same section in report.
    Only works when `each` is not specified.
- `cases` (`type=json`): If you have multiple cases, you can specify them here.
    Keys are the names of the cases and values are the options above except `mutaters`.
    If some options are not specified, the options in `envs` will be used.
    If no cases are specified, a default case will be used with case name `DEFAULT`.

## Examples

```toml
[CellsDistribution.envs.mutaters]
# Add Patient1_Tumor_Expanded column with CDR3.aa that
# expands in Tumor of patient 1
Patient1_Tumor_Expanded = '''
  expanded(., region, "Tumor", subset = patient == "Lung1", uniq = FALSE)
'''

[CellsDistribution.envs.cases.Patient1_Tumor_Expanded]
cells_by = "Patient1_Tumor_Expanded"
cells_orderby = "desc(CloneSize)"
group_by = "region"
group_order = [ "Tumor", "Normal" ]
```

![CellsDistribution_example](images/CellsDistribution_example.png)
