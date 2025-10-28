# MetabolicFeatures

This process performs enrichment analysis for the metabolic pathways for each group in each subset.

The enrichment analysis is done with [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
package or the [`GSEA_R`](https://github.com/GSEA-MSigDB/GSEA_R) package.<br />

## Input

- `sobjfile`:
    The Seurat object file in rds.<br />
    It should be loaded as a Seurat object

## Output

- `outdir`: *Default: `{{in.sobjfile | stem}}.pathwayfeatures`*. <br />
    The output directory.<br />
    It will contain the GSEA results and plots.<br />

## Environment Variables

- `ncores` *(`type=int;pgarg`)*: *Default: `1`*. <br />
    Number of cores to use for parallelization for
    the comparisons for each subset and group.<br />
    Defaults to `ScrnaMetabolicLandscape.ncores`.<br />
- `prerank_method` *(`choice`)*: *Default: `signal_to_noise`*. <br />
    Method to use for gene preranking.<br />
    Signal to noise: the larger the differences of the means
    (scaled by the standard deviations); that is, the more distinct
    the gene expression is in each phenotype and the more the gene
    acts as a “class marker.”.<br />
    Absolute signal to noise: the absolute value of the signal to
    noise.<br />
    T test: Uses the difference of means scaled by the standard
    deviation and number of samples.<br />
    Ratio of classes: Uses the ratio of class means to calculate
    fold change for natural scale data.<br />
    Diff of classes: Uses the difference of class means to calculate
    fold change for nature scale data
    Log2 ratio of classes: Uses the log2 ratio of class means to
    calculate fold change for natural scale data. This is the
    recommended statistic for calculating fold change for log scale
    data.<br />
    - `signal_to_noise`:
        Signal to noise
    - `s2n`:
        Alias of signal_to_noise
    - `abs_signal_to_noise`:
        absolute signal to noise
    - `abs_s2n`:
        Alias of abs_signal_to_noise
    - `t_test`:
        T test
    - `ratio_of_classes`:
        Also referred to as fold change
    - `diff_of_classes`:
        Difference of class means
    - `log2_ratio_of_classes`:
        Log2 ratio of class means
- `gmtfile` *(`pgarg`)*:
    The GMT file with the metabolic pathways.<br />
    Defaults to `ScrnaMetabolicLandscape.gmtfile`
- `subset_by` *(`pgarg;readonly`)*:
    Subset the data by the given column in the
    metadata. For example, `Response`.<br />
    `NA` values will be removed in this column.<br />
    Defaults to `ScrnaMetabolicLandscape.subset_by`
    If None, the data will not be subsetted.<br />
- `group_by` *(`pgarg;readonly`)*:
    Group the data by the given column in the
    metadata. For example, `cluster`.<br />
    Defaults to `ScrnaMetabolicLandscape.group_by`
- `comparisons` *(`type=list`)*: *Default: `[]`*. <br />
    The comparison groups to use for the analysis.<br />
    If not provided, each group in the `group_by` column will be used
    to compare with the other groups.<br />
    If a single group is provided as an element, it will be used to
    compare with all the other groups.<br />
    For example, if we have `group_by = "cluster"` and we have
    `1`, `2` and `3` in the `group_by` column, we could have
    `comparisons = ["1", "2"]`, which will compare the group `1` with groups
    `2` and `3`, and the group `2` with groups `1` and `3`. We could also
    have `comparisons = ["1:2", "1:3"]`, which will compare the group `1` with
    group `2` and group `1` with group `3`.<br />
- `fgsea_args` *(`type=json`)*: *Default: `{}`*. <br />
    Other arguments for the `fgsea::fgsea()` function.<br />
    For example, `{"minSize": 15, "maxSize": 500}`.<br />
    See <https://rdrr.io/bioc/fgsea/man/fgsea.html> for more details.<br />
- `plots` *(`type=json`)*: *Default: `{'Summary Plot': Diot({'plot_type': 'summary', 'top_term': 10, 'devpars': Diot({'res': 100})}), 'Enrichment Plots': Diot({'plot_type': 'gsea', 'top_term': 10, 'devpars': Diot({'res': 100})})}`*. <br />
    The plots to generate.<br />
    Names will be used as the title for the plot. Values will be the arguments
    passed to `biopipen.utils::VizGSEA()` function.<br />
    See <https://user.github.io/biopipen.utils.R/reference/VizGSEA.html>.<br />
    A key `level` is supported to specify the level of the plot.<br />
    Possible values are `case`, which includes all subsets and groups in the
    case; `subset`, which includes all groups in the subset; otherwise, it
    will plot for the groups.<br />
    For `case`/`subset` level plots, current `plot_type` only "dot" is supported
    for now, then the values will be passed to `plotthis::DotPlot()`
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    Multiple cases for the analysis.<br />
    If you only have one case, you can specify the parameters directly to
    `envs.prerank_method`, `envs.subset_by`, `envs.group_by`,
    `envs.comparisons`, `envs.fgsea_args` and `envs.plots`.<br />
    The name of this default case will be `envs.subset_by`.<br />
    If you have multiple cases, you can specify the parameters for each case
    in a dictionary. The keys will be the names of the cases and the values
    will be dictionaries with the parameters for each case, where the values
    will be inherited from `envs.prerank_method`,
    `envs.subset_by`, `envs.group_by`, `envs.comparisons`, `envs.fgsea_args`
    and `envs.plots`.<br />

