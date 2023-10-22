# ScFGSEA

Gene set enrichment analysis for cells in different groups using [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html).

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
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`, that can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).<br />
    For example, you can use `{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1')"}`
    to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
    with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1. The values in this columns for other clones will be `NA`.<br />
    Those functions take following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `group-by`: The column name in metadata to group the cells.<br />
    * `idents`: The first group or both groups of cells to compare (value in `group-by` column). If only the first group is given, the rest of the cells (with non-NA in `group-by` column) will be used as the second group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`)
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `order`: The order of the returned ids. It could be `sum` or `diff`, which is the sum or diff of the `compare` between idents.<br />
    Two kinds of modifiers can be added, including `desc` and `abs`.<br />
    For example, `sum,desc` means the sum of `compare` between idents in descending order.<br />
    Default is `diff,abs,desc`. It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
    ids will be in the same order as in `df`.<br />
    Note that the numeric column should be the same for all cells in the same group. This will not be checked (only the first value is used).<br />
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
    The pathways in GMT format, with the gene names/ids in the same format as the seurat object
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

