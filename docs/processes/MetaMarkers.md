# MetaMarkers

Find markers between three or more groups of cells, using one-way ANOVA or Kruskal-Wallis test.

Sometimes, you may want to find the markers for cells from more than 2 groups.<br />
In this case, you can use this process to find the markers for the groups and
do enrichment analysis for the markers. Each marker is examined using either
one-way ANOVA or Kruskal-Wallis test.<br />
The p values are adjusted using the specified method. The significant markers
are then used for enrichment analysis using
[enrichr](https://maayanlab.cloud/Enrichr/) api.<br />

Other than the markers and the enrichment analysis as outputs, this process also
generates violin plots for the top 10 markers.<br />

## Environment Variables

- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    Number of cores to use to parallelize for genes
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata
    The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.<br />
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`,
    which can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).<br />
    For example, you can use
    `{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq = FALSE)"}`
    to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
    with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1.<br />
    The values in this columns for other clones will be `NA`.<br />
    Those functions take following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `group-by`: The column name in metadata to group the cells.<br />
    * `idents`: The first group or both groups of cells to compare (value in `group-by` column). If only the first group is given, the rest of the cells (with non-NA in `group-by` column) will be used as the second group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).<br />
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    If numeric column is given, the values should be the same for all cells in the same group.<br />
    This will not be checked (only the first value is used).<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `order`: The order of the returned ids. It could be `sum` or `diff`, which is the sum or diff of the `compare` between idents.<br />
    Two kinds of modifiers can be added, including `desc` and `abs`.<br />
    For example, `sum,desc` means the sum of `compare` between idents in descending order.<br />
    Default is `diff,abs,desc`. It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
    ids will be in the same order as in `df`.<br />
    * `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.<br />
    * `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.<br />
- `group-by`:
    The column name in metadata to group the cells.<br />
    If only `group-by` is specified, and `idents` are
    not specified, markers will be found for all groups in this column.<br />
    `NA` group will be ignored.<br />
- `idents`:
    The groups of cells to compare, values should be in the `group-by` column.<br />
- `each`:
    The column name in metadata to separate the cells into different cases.<br />
- `prefix_each` *(`flag`)*: *Default: `True`*. <br />
    Whether to add the `each` value as prefix to the case name.<br />
- `dbs` *(`list`)*: *Default: `['GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021', 'KEGG_2021_Human']`*. <br />
    The dbs to do enrichment analysis for significant
    markers See below for all libraries.<br />
    <https://maayanlab.cloud/Enrichr/#libraries>
- `subset`:
    The subset of the cells to do the analysis.<br />
    An expression passed to `dplyr::filter()`.<br />
- `p_adjust` *(`choice`)*: *Default: `BH`*. <br />
    The method to adjust the p values, which can be used to filter the significant markers.<br />
    See also <https://rdrr.io/r/stats/p.adjust.html>
    - `holm`:
        Holm-Bonferroni method
    - `hochberg`:
        Hochberg method
    - `hommel`:
        Hommel method
    - `bonferroni`:
        Bonferroni method
    - `BH`:
        Benjamini-Hochberg method
    - `BY`:
        Benjamini-Yekutieli method
    - `fdr`:
        FDR method of Benjamini-Hochberg
    - `none`:
        No adjustment
- `sigmarkers`: *Default: `p_adjust < 0.05`*. <br />
    An expression passed to `dplyr::filter()` to filter the
    significant markers for enrichment analysis. The default is `p.value < 0.05`.<br />
    If `method = 'anova'`, the variables that can be used for filtering are:<br />
    `sumsq`, `meansq`, `statistic`, `p.value` and `p_adjust`.<br />
    If `method = 'kruskal'`, the variables that can be used for filtering are:<br />
    `statistic`, `p.value` and `p_adjust`.<br />
- `section`: *Default: `DEFAULT`*. <br />
    The section name for the report.<br />
    Worked only when `each` is not specified.<br />
    Otherwise, the section name will be constructed from `each` and `group-by`.<br />
    If `DEFAULT`, and it's the only section, it not included in the case/section names.<br />
- `method` *(`choice`)*: *Default: `anova`*. <br />
    The method for the test.<br />
    - `anova`:
        One-way ANOVA
    - `kruskal`:
        Kruskal-Wallis test
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    If you have multiple cases, you can specify them
    here. The keys are the names of the cases and the values are the
    above options except `ncores` and `mutaters`. If some options are
    not specified, the default values specified above will be used.<br />
    If no cases are specified, the default case will be added with
    the default values under `envs` with the name `DEFAULT`.<br />

