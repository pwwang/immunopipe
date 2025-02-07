# CellsDistribution

Distribution of cells (i.e. in a TCR clone) from different groups for each cluster

This generates a set of pie charts with proportion of cells in each cluster
Rows are the cells identities (i.e. TCR clones or TCR clusters), columns
are groups (i.e. clinic groups).<br />

## Input

- `srtobj`:
    The seurat object in RDS format

## Output

- `outdir`: *Default: `{{in.srtobj | stem}}.cells_distribution`*. <br />
    The output directory

## Environment Variables

- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata
    Keys are the names of the mutaters and values are the R expressions
    passed by `dplyr::mutate()` to mutate the metadata.<br />
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`,
    which can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).<br />
    See also <../configurations/#mutater-helpers>.<br />
    For example, you can use
    `{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq = FALSE)"}`
    to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
    with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1.<br />
    The values in this columns for other clones will be `NA`.<br />
    Those functions take following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `group.by`: The column name in metadata to group the cells.<br />
    * `idents`: The first group or both groups of cells to compare (value in `group.by` column). If only the first group is given, the rest of the cells (with non-NA in `group.by` column) will be used as the second group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `each`: A column name (without quotes) in metadata to split the cells.<br />
    Each comparison will be done for each value in this column (typically each patient or subject).<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).<br />
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    If numeric column is given, the values should be the same for all cells in the same group.<br />
    This will not be checked (only the first value is used).<br />
    It is helpful to use `Clones` to use the raw clone size from TCR data, in case the cells are not completely mapped to RNA data.<br />
    Also if you have `subset` set or `NA`s in `group.by` column, you should use `.n` to compare the number of cells in each group.<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.<br />
    * `order`: The expression passed to `dplyr::arrange()` to order intermediate dataframe and get the ids in order accordingly.<br />
    The intermediate dataframe includes the following columns:<br />
    * `<id>`: The ids of clones (i.e. `CDR3.aa`).<br />
    * `<each>`: The values in `each` column.<br />
    * `ident_1`: The size of clones in the first group.<br />
    * `ident_2`: The size of clones in the second group.<br />
    * `.diff`: The difference between the sizes of clones in the first and second groups.<br />
    * `.sum`: The sum of the sizes of clones in the first and second groups.<br />
    * `.predicate`: Showing whether the clone is expanded/collapsed/emerged/vanished.<br />
    * `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.<br />
    * `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.<br />

    You can also use `top()` to get the top clones (i.e. the clones with the largest size) in each group.<br />
    For example, you can use
    `{"Patient1_Top10_Clones": "top(subset = Patent == 'Patient1', uniq = FALSE)"}`
    to create a new column in metadata named `Patient1_Top10_Clones`.<br />
    The values in this columns for other clones will be `NA`.<br />
    This function takes following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).<br />
    * `n`: The number of top clones to return. Default is `10`.<br />
    If n < 1, it will be treated as the percentage of the size of the group.<br />
    Specify `0` to get all clones.<br />
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    If numeric column is given, the values should be the same for all cells in the same group.<br />
    This will not be checked (only the first value is used).<br />
    It is helpful to use `Clones` to use the raw clone size from TCR data, in case the cells are not completely mapped to RNA data.<br />
    Also if you have `subset` set or `NA`s in `group.by` column, you should use `.n` to compare the number of cells in each group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `each`: A column name (without quotes) in metadata to split the cells.<br />
    Each comparison will be done for each value in this column (typically each patient or subject).<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.<br />
    * `with_ties`: Whether to include ties (i.e. clones with the same size as the last clone) or not. Default is `FALSE`.<br />

- `cluster_orderby`:
    The order of the clusters to show on the plot.<br />
    An expression passed to `dplyr::summarise()` on the grouped data frame (by `seurat_clusters`).<br />
    The summary stat will be passed to `dplyr::arrange()` to order the clusters. It's applied on the whole meta.data before grouping and subsetting.<br />
    For example, you can order the clusters by the activation score of
    the cluster: `desc(mean(ActivationScore, na.rm = TRUE))`, suppose you have a column
    `ActivationScore` in the metadata.<br />
- `group_by`:
    The column name in metadata to group the cells for the columns of the plot.<br />
- `group_order` *(`list`)*: *Default: `[]`*. <br />
    The order of the groups (columns) to show on the plot
- `cells_by`:
    The column name in metadata to group the cells for the rows of the plot.<br />
    If your cell groups have overlapping cells, you can also use multiple columns, separated by comma (`,`).<br />
    These columns will be concatenated to form the cell groups. For the overlapping cells, they will be
    counted multiple times for different groups. So make sure the cell group names in different columns
    are unique.<br />
- `cells_order` *(`list`)*: *Default: `[]`*. <br />
    The order of the cells (rows) to show on the plot
- `cells_orderby`:
    An expression passed to `dplyr::arrange()` to order the cells (rows) of the plot.<br />
    Only works when `cells-order` is not specified.<br />
    The data frame passed to `dplyr::arrange()` is grouped by `cells_by` before ordering.<br />
    You can have multiple expressions separated by semicolon (`;`). The expessions will be parsed by `rlang::parse_exprs()`.<br />
    4 extra columns were added to the metadata for ordering the rows in the plot:<br />
    * `CloneSize`: The size (number of cells) of clones (identified by `cells_by`)
    * `CloneGroupSize`: The clone size in each group (identified by `group_by`)
    * `CloneClusterSize`: The clone size in each cluster (identified by `seurat_clusters`)
    * `CloneGroupClusterSize`: The clone size in each group and cluster (identified by `group_by` and `seurat_clusters`)
- `cells_n` *(`type=int`)*: *Default: `10`*. <br />
    The max number of groups to show for each cell group identity (row).<br />
    Ignored if `cells_order` is specified.<br />
- `subset`:
    An expression to subset the cells, will be passed to `dplyr::filter()` on metadata.<br />
    This will be applied prior to `each`.<br />
- `descr`:
    The description of the case, will be shown in the report.<br />
- `hm_devpars` *(`ns`)*:
    The device parameters for the heatmaps.<br />
    - `res` *(`type=int`)*:
        The resolution of the heatmaps.<br />
    - `height` *(`type=int`)*:
        The height of the heatmaps.<br />
    - `width` *(`type=int`)*:
        The width of the heatmaps.<br />
- `devpars` *(`ns`)*:
    The device parameters for the plots of pie charts.<br />
    - `res` *(`type=int`)*:
        The resolution of the plots
    - `height` *(`type=int`)*:
        The height of the plots
    - `width` *(`type=int`)*:
        The width of the plots
- `each`:
    The column name in metadata to separate the cells into different plots.<br />
- `prefix_each` *(`flag`)*: *Default: `True`*. <br />
    Whether to prefix the `each` column name to the
    value as the case/section name.<br />
- `section`: *Default: `DEFAULT`*. <br />
    The section to show in the report. This allows different cases to be put in the same section in report.<br />
    Only works when `each` is not specified.<br />
    The `section` is used to collect cases and put the results under the same directory and the same section in report.<br />
    When `each` for a case is specified, the `section` will be ignored and case name will be used as `section`.<br />
    The cases will be the expanded values in `each` column. When `prefix_each` is True, the column name specified by `each` will be prefixed to each value as directory name and expanded case name.<br />
- `overlap` *(`list`)*: *Default: `[]`*. <br />
    Plot the overlap of cell groups (values of `cells_by`) in different cases
    under the same section.<br />
    The section must have at least 2 cases, each case should have a single `cells_by` column.<br />
- `cases` *(`type=json;order=99`)*: *Default: `{}`*. <br />
    If you have multiple cases, you can specify them here.<br />
    Keys are the names of the cases and values are the options above except `mutaters`.<br />
    If some options are not specified, the options in `envs` will be used.<br />
    If no cases are specified, a default case will be used with case name `DEFAULT`.<br />

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

![CellsDistribution_example](../latest/processes/images/CellsDistribution_example.png)

