# MetabolicFeatures

This process performs enrichment analysis for the metabolic pathways for each group in each subset.

The enrichment analysis is done with [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
package or the [`GSEA_R`](https://github.com/GSEA-MSigDB/GSEA_R) package.<br />

## Input

- `sobjfile`:

## Output

- `outdir`: *Default: `{{in.sobjfile | stem}}.pathwayfeatures`*. <br />

## Environment Variables

- `ncores` *(`type=int;pgarg`)*: *Default: `1`*. <br />
    Number of cores to use for parallelization.<br />
    Defaults to `ScrnaMetabolicLandscape.ncores`
- `fgsea` *(`flag`)*: *Default: `True`*. <br />
    Whether to do fast gsea analysis using `fgsea` package.<br />
    If `False`, the `GSEA_R` package will be used.<br />
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
- `top` *(`type=int`)*: *Default: `10`*. <br />
    N top of enriched pathways to show
- `gmtfile` *(`pgarg`)*:
    The GMT file with the metabolic pathways.<br />
    Defaults to `ScrnaMetabolicLandscape.gmtfile`
- `grouping` *(`type=auto;pgarg;readonly`)*:
    Defines the basic groups to
    investigate the metabolic activity.<br />
    Defaults to `ScrnaMetabolicLandscape.grouping`
- `grouping_prefix` *(`type=auto;pgarg;readonly`)*: *Default: `""`*. <br />
    Working as a prefix to
    group names.<br />
    Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
- `subsetting` *(`type=auto;pgarg;readonly`)*:
    How do we subset the data.<br />
    Another column(s) in the metadata.<br />
    Defaults to `ScrnaMetabolicLandscape.subsetting`
- `subsetting_prefix` *(`type=auto;pgarg;readonly`)*:
    Working as a prefix to
    subset names.<br />
    Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

