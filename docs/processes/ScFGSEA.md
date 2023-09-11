# ScFGSEA

This process allows us to do Gene Set Enrichment Analysis (GSEA) on the expression data, but based on variaties of grouping, including the from the meta data and the scTCR-seq data as well.

The GSEA is done using the [fgsea][1] package, which allows to quickly and accurately calculate arbitrarily low GSEA P-values for a collection of gene sets. The fgsea package is based on the fast algorithm for preranked GSEA described in [Subramanian et al. 2005](https://www.pnas.org/content/102/43/15545).

For each case, the process will generate a table with the enrichment scores for each gene set, and GSEA plots for the top gene sets.

## Environment variables

- `ncores` (`type=int`): Number of cores for parallelization
    Passed to `nproc` of [`fgseaMultilevel()`][2].
- `mutaters` (`type=json`): The mutaters to mutate the metadata.
    The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`, that can be used to identify the expanded/collapsed/emerged/vanished groups (i.e. TCR clones).
    For example, you can use `{"Patient1_Tumor_Expanded_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq=FALSE)"}`
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
- `group-by`: The column name in metadata to group the cells.
- `ident-1`: The first group of cells to compare
- `ident-2`: The second group of cells to compare, if not provided, the rest of the cells that are not `NA`s in `group-by` column are used for `ident-2`.
- `each`: The column name in metadata to separate the cells into different subsets to do the analysis.
- `section`: The section name for the report. Worked only when `each` is not specified. Otherwise, the section name will be constructed from `each` and its value.
    This allows different cases to be put into the same section in the report.
- `gmtfile`: The pathways in GMT format, with the gene names/ids in the same format as the seurat object
- `method` (`choice`): The method to do the preranking.
    - `signal_to_noise`: Signal to noise.
        The larger the differences of the means (scaled by the standard deviations);
        that is, the more distinct the gene expression is in each phenotype and the more the gene
        acts as a "class marker".
    - `s2n`: Alias of signal_to_noise.
    - `abs_signal_to_noise`: The absolute value of signal_to_noise.
    - `abs_s2n`: Alias of abs_signal_to_noise.
    - `t_test`: T test.
        Uses the difference of means scaled by the standard deviation and number of samples.
    - `ratio_of_classes`: Also referred to as fold change.
        Uses the ratio of class means to calculate fold change for natural scale data.
    - `diff_of_classes`: Difference of class means.
        Uses the difference of class means to calculate fold change for nature scale data
    - `log2_ratio_of_classes`: Log2 ratio of class means.
        Uses the log2 ratio of class means to calculate fold change for natural scale data.
        This is the recommended statistic for calculating fold change for log scale data.
- `top` (`type=auto`): Do gsea table and enrich plot for top N pathways.
    If it is < 1, will apply it to `padj`, selecting pathways with `padj` < `top`.
- `eps` (`type=float`): This parameter sets the boundary for calculating the p value.
    See <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
- `minsize` (`type=int`): Minimal size of a gene set to test. All pathways below the threshold are excluded.
- `maxsize` (`type=int`): Maximal size of a gene set to test. All pathways above the threshold are excluded.
- `rest` (`type=json`): Rest arguments for [`fgsea()`](https://rdrr.io/bioc/fgsea/man/fgsea.html)
    See also <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
- `cases` (`type=json`): If you have multiple cases, you can specify them here.
    The keys are the names of the cases and the values are the above options except `mutaters`.
    If some options are not specified, the default values specified above will be used.
    If no cases are specified, the default case will be added with the name `DEFAULT`.

[1]: https://bioconductor.org/packages/release/bioc/html/fgsea.html
[2]: https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html
