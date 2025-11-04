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

![MetabolicPathwayHeterogeneity](../latest/processes/images/MetabolicPathwayHeterogeneity.png)

## Input

- `sobjfile`:

## Output

- `outdir`: *Default: `{{in.sobjfile | stem}}.pathwayhetero`*. <br />

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
- `fgsea_args` *(`type=json`)*: *Default: `{'scoreType': 'std', 'nproc': 1}`*. <br />
    Other arguments for the `fgsea::fgsea()` function.<br />
    For example, `{"minSize": 15, "maxSize": 500}`.<br />
    See <https://rdrr.io/bioc/fgsea/man/fgsea.html> for more details.<br />
- `plots` *(`type=json`)*: *Default: `{'Pathway Heterogeneity': Diot({'plot_type': 'dot', 'devpars': Diot({'res': 100})})}`*. <br />
    The plots to generate.<br />
    Names will be used as the title for the plot. Values will be the arguments
    passed to `biopipen.utils::VizGSEA()` function.<br />
    See <https://pwwang.github.io/biopipen.utils.R/reference/VizGSEA.html>.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    Multiple cases for the analysis.<br />
    If you only have one case, you can specify the parameters directly to
    `envs.subset_by`, `envs.group_by`, `envs.fgsea_args`, `envs.plots`,
    `envs.select_pcs`, and `envs.pathway_pval_cutoff`.<br />
    The name of this default case will be `envs.subset_by`.<br />
    If you have multiple cases, you can specify the parameters for each case
    in a dictionary. The keys will be the names of the cases and the values
    will be dictionaries with the parameters for each case, where the values
    will be inherited from `envs.subset_by`, `envs.group_by`, `envs.fgsea_args`,
    `envs.plots`, `envs.select_pcs`, and `envs.pathway_pval_cutoff`.<br />

