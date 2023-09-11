# MetabolicFeatures

This process performs enrichment analysis for the metabolic pathways for each group in each subset. The enrichment analysis is done with [`fgsea`][1] package or the [`GSEA_R`][2] package.

## Environment variables

- `ncores` (`type=int`): Number of cores to use for parallelization
    Defaults to `ScrnaMetabolicLandscape.ncores`
- `fgsea` (`flag`): Whether to do fast gsea analysis using [`fgsea`][1] package.
    If `False`, the [`GSEA_R`][2] package will be used.
- `prerank_method` (`choice`): Method to use for gene preranking if [`fgsea`][1] is used.
    Signal to noise: the larger the differences of the means
    (scaled by the standard deviations); that is, the more distinct
    the gene expression is in each phenotype and the more the gene
    acts as a “class marker.”.
    Absolute signal to noise: the absolute value of the signal to
    noise.
    T test: Uses the difference of means scaled by the standard
    deviation and number of samples.
    Ratio of classes: Uses the ratio of class means to calculate
    fold change for natural scale data.
    Diff of classes: Uses the difference of class means to calculate
    fold change for nature scale data
    Log2 ratio of classes: Uses the log2 ratio of class means to
    calculate fold change for natural scale data. This is the
    recommended statistic for calculating fold change for log scale
    data.
    - `signal_to_noise`: Signal to noise
    - `s2n`: Alias of signal_to_noise
    - `abs_signal_to_noise`: absolute signal to noise
    - `abs_s2n`: Alias of abs_signal_to_noise
    - `t_test`: T test
    - `ratio_of_classes`: Also referred to as fold change
    - `diff_of_classes`: Difference of class means
    - `log2_ratio_of_classes`: Log2 ratio of class means
- `top` (`type=int`): Top N of enriched pathways to show
- `gmtfile`: The GMT file with the metabolic pathways.
    Defaults to `ScrnaMetabolicLandscape.gmtfile`
- `grouping` (`type=auto`;`readonly`): Defines the basic groups to
    investigate the metabolic activity.
    Defaults to `ScrnaMetabolicLandscape.grouping`
- `grouping_prefix` (`type=auto`;`readonly`): Working as a prefix to
    group names.
    Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
- `subsetting` (`type=auto`;`readonly`): How do we subset the data.
    Another column(s) in the metadata.
    Defaults to `ScrnaMetabolicLandscape.subsetting`
- `subsetting_prefix` (`type=auto`;`readonly`): Working as a prefix to
    subset names.
    Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

[1]: https://bioconductor.org/packages/release/bioc/html/fgsea.html
[2]: https://github.com/GSEA-MSigDB/GSEA_R
