# ScFGSEA

Gene set enrichment analysis for cells in different groups using `fgsea`

This process allows us to do Gene Set Enrichment Analysis (GSEA) on the expression data,
but based on variaties of grouping, including the from the meta data and the
scTCR-seq data as well.<br />

The GSEA is done using the
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) package,
which allows to quickly and accurately calculate arbitrarily low GSEA P-values
for a collection of gene sets.<br />
The fgsea package is based on the fast algorithm for preranked GSEA described in
[Subramanian et al. 2005](https://www.pnas.org/content/102/43/15545).<br />

For each case, the process will generate a table with the enrichment scores for
each gene set, and GSEA plots for the top gene sets.<br />

## Environment Variables

- `ncores` *(`type=int`)*: *Default: `1`*. <br />
    Number of cores for parallelization
    Passed to `nproc` of `fgseaMultilevel()`.<br />
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    The mutaters to mutate the metadata.<br />
    The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.<br />
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
- `group-by`:
    The column name in metadata to group the cells.<br />
- `ident-1`:
    The first group of cells to compare
- `ident-2`:
    The second group of cells to compare, if not provided, the rest of the cells that are not `NA`s in `group-by` column are used for `ident-2`.<br />
- `each`:
    The column name in metadata to separate the cells into different subsets to do the analysis.<br />
- `section`: *Default: `DEFAULT`*. <br />
    The section name for the report. Worked only when `each` is not specified. Otherwise, the section name will be constructed from `each` and its value.<br />
    This allows different cases to be put into the same section in the report.<br />
- `gmtfile`: *Default: `""`*. <br />
    The pathways in GMT format, with the gene names/ids in the same format as the seurat object.<br />
    One could also use a URL to a GMT file. For example, from <https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Pathways/>.<br />
- `method` *(`choice`)*: *Default: `s2n`*. <br />
    The method to do the preranking.<br />
    - `signal_to_noise`:
        Signal to noise.<br />
        The larger the differences of the means (scaled by the standard deviations);
        that is, the more distinct the gene expression is in each phenotype and the more the gene
        acts as a "class marker".<br />
    - `s2n`:
        Alias of signal_to_noise.<br />
    - `abs_signal_to_noise`:
        The absolute value of signal_to_noise.<br />
    - `abs_s2n`:
        Alias of abs_signal_to_noise.<br />
    - `t_test`:
        T test.<br />
        Uses the difference of means scaled by the standard deviation and number of samples.<br />
    - `ratio_of_classes`:
        Also referred to as fold change.<br />
        Uses the ratio of class means to calculate fold change for natural scale data.<br />
    - `diff_of_classes`:
        Difference of class means.<br />
        Uses the difference of class means to calculate fold change for nature scale data
    - `log2_ratio_of_classes`:
        Log2 ratio of class means.<br />
        Uses the log2 ratio of class means to calculate fold change for natural scale data.<br />
        This is the recommended statistic for calculating fold change for log scale data.<br />
- `top` *(`type=auto`)*: *Default: `20`*. <br />
    Do gsea table and enrich plot for top N pathways.<br />
    If it is < 1, will apply it to `padj`, selecting pathways with `padj` < `top`.<br />
- `eps` *(`type=float`)*: *Default: `0`*. <br />
    This parameter sets the boundary for calculating the p value.<br />
    See <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
- `minsize` *(`type=int`)*: *Default: `10`*. <br />
    Minimal size of a gene set to test. All pathways below the threshold are excluded.<br />
- `maxsize` *(`type=int`)*: *Default: `100`*. <br />
    Maximal size of a gene set to test. All pathways above the threshold are excluded.<br />
- `rest` *(`type=json;order=98`)*: *Default: `{}`*. <br />
    Rest arguments for [`fgsea()`](https://rdrr.io/bioc/fgsea/man/fgsea.html)
    See also <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
- `cases` *(`type=json;order=99`)*: *Default: `{}`*. <br />
    If you have multiple cases, you can specify them here.<br />
    The keys are the names of the cases and the values are the above options except `mutaters`.<br />
    If some options are not specified, the default values specified above will be used.<br />
    If no cases are specified, the default case will be added with the name `DEFAULT`.<br />

