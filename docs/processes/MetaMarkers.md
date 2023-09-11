# MetaMarkers

Sometimes, you may want to find the markers for cells from more than 2 groups. In this case, you can use this process to find the markers for the groups and do enrichment analysis for the markers. Each marker is examined using either one-way ANOVA or Kruskal-Wallis test. The p values are adjusted using the specified method. The significant markers are then used for enrichment analysis using [enrichr](https://maayanlab.cloud/Enrichr/) api.

Other than the markers and the enrichment analysis as outputs, this process also generates violin plots for the top 10 markers.

## Environment variables

- `ncores` (`type=int`): Number of cores to use to parallelize for genes
- `mutaters` (`type=json`): The mutaters to mutate the metadata
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`, that can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).
    For example, you can use `{"Patient1_Tumor_Collapsed_Clones": "expanded(Source, 'Tumor', subset = Patent == 'Patient1')"}`
    to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
    with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1. The values in this columns for other clones will be `NA`.
    Those functions take following arguments:
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
`group-by`: The column name in metadata to group the cells.
    If only `group-by` is specified, and `idents` are
    not specified, markers will be found for all groups in this column.
    `NA` group will be ignored.
- `idents`: The groups of cells to compare, values should be in the `group-by` column.
- `each`: The column name in metadata to separate the cells into different cases.
- `prefix_each` (`flag`): Whether to add the `each` value as prefix to the case name.
- `dbs` (`list`): The dbs to do enrichment analysis for significant
    markers See below for all libraries.
    <https://maayanlab.cloud/Enrichr/#libraries>
- `p_adjust` (`choice`): The method to adjust the p values, which can be used to filter the significant markers.
    See also <https://rdrr.io/r/stats/p.adjust.html>
    - `holm`: Holm-Bonferroni method
    - `hochberg`: Hochberg method
    - `hommel`: Hommel method
    - `bonferroni`: Bonferroni method
    - `BH`: Benjamini-Hochberg method
    - `BY`: Benjamini-Yekutieli method
    - `fdr`: FDR method of Benjamini-Hochberg
    - `none`: No adjustment
- `sigmarkers`: An expression passed to `dplyr::filter()` to filter the
    significant markers for enrichment analysis. The default is `p.value < 0.05`.
    If `method = 'anova'`, the variables that can be used for filtering are:
    `sumsq`, `meansq`, `statistic`, `p.value` and `p_adjust`.
    If `method = 'kruskal'`, the variables that can be used for filtering are:
    `statistic`, `p.value` and `p_adjust`.
- `section`: The section name for the report.
    Worked only when `each` is not specified.
    Otherwise, the section name will be constructed from `each` and `group-by`.
    If `DEFAULT`, and it's the only section, it not included in the case/section names.
- `method` (`choice`): The method for the test.
    - `anova`: One-way ANOVA
    - `kruskal`: Kruskal-Wallis test
- `cases` (`type=json`): If you have multiple cases, you can specify them
    here. The keys are the names of the cases and the values are the
    above options except `ncores` and `mutaters`. If some options are
    not specified, the default values specified above will be used.
    If no cases are specified, the default case will be added with
    the default values under `envs` with the name `DEFAULT`.
