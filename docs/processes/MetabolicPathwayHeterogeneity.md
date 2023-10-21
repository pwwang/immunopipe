# MetabolicPathwayHeterogeneity

Calculate Metabolic Pathway heterogeneity.

For each subset, the normalized enrichment score (NES) of each metabolic pathway
is calculated for each group.<br />
The NES is calculated by comparing the enrichment score of the subset to the
enrichment scores of the same subset in the permutations.<br />
The p-value is calculated by comparing the NES to the NESs of the same subset
in the permutations.<br />
The heterogeneity can be reflected by the NES values and the p-values in
different groups for the metabolic pathways.<br />

![MetabolicPathwayHeterogeneity](https://pwwang.github.io/immunopipe/processes/images/MetabolicPathwayHeterogeneity.png)

## Environment Variables

- `gmtfile` *(`pgarg`)*:
    The GMT file with the metabolic pathways.<br />
    Defaults to `ScrnaMetabolicLandscape.gmtfile`
- `select_pcs` *(`type=float`)*: *Default: `0.8`*. <br />
    Select the PCs to use for the analysis.<br />
- `pathway_pval_cutoff` *(`type=float`)*: *Default: `0.01`*. <br />
    The p-value cutoff to select
    the enriched pathways
- `ncores` *(`type=int;pgarg`)*: *Default: `1`*. <br />
    Number of cores to use for parallelization
    Defaults to `ScrnaMetabolicLandscape.ncores`
- `bubble_devpars` *(`ns`)*:
    The devpars for the bubble plot
    - `width` *(`type=int`)*:
        The width of the plot
    - `height` *(`type=int`)*:
        The height of the plot
    - `res` *(`type=int`)*:
        The resolution of the plot
- `grouping` *(`type=auto;pgarg;readonly`)*:
    Defines the basic groups to
    investigate the metabolic activity.<br />
    Defaults to `ScrnaMetabolicLandscape.grouping`
- `grouping_prefix` *(`type=auto;pgarg;readonly`)*: *Default: `""`*. <br />
    Working as a prefix to group
    names.<br />
    Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
- `subsetting` *(`type=auto;pgarg;readonly`)*:
    How do we subset the data.<br />
    Another column(s) in the metadata.<br />
    Defaults to `ScrnaMetabolicLandscape.subsetting`
- `subsetting_prefix` *(`type=auto;pgarg;readonly`)*:
    Working as a prefix to
    subset names.<br />
    Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

