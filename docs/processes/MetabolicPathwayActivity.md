# MetabolicPathwayActivity

This process calculates the pathway activities in different groups and subsets.

The cells are first grouped by subsets and then the metabolic activities are
examined for each groups in different subsets.<br />

For each subset, a heatmap and a violin plot will be generated.<br />
The heatmap shows the pathway activities for each group and each metabolic pathway

![MetabolicPathwayActivity_heatmap](images/MetabolicPathwayActivity_heatmap.png){: width="80%"}

The violin plot shows the distribution of the pathway activities for each group

![MetabolicPathwayActivity_violin](images/MetabolicPathwayActivity_violin.png){: width="45%"}

## Input

- `sobjfile`:
    The Seurat object file.<br />
    It should be loaded as a Seurat object

## Output

- `outdir`: *Default: `{{in.sobjfile | stem}}.pathwayactivity`*. <br />
    The output directory.<br />
    It will contain the pathway activity score files and plots.<br />

## Environment Variables

- `ntimes` *(`type=int`)*: *Default: `5000`*. <br />
    Number of permutations to estimate the p-values
- `ncores` *(`type=int;pgarg`)*: *Default: `1`*. <br />
    Number of cores to use for parallelization
    Defaults to `ScrnaMetabolicLandscape.ncores`
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
- `plots` *(`type=json`)*: *Default: `{'Pathway Activity (violin plot)': Diot({'plot_type': 'violin', 'add_box': True, 'devpars': Diot({'res': 100})}), 'Pathway Activity (heatmap)': Diot({'plot_type': 'heatmap', 'devpars': Diot({'res': 100})})}`*. <br />
    The plots to generate.<br />
    Names will be used as the prefix for the output files. Values will be
    a dictionary with the following keys:<br />
    * `plot_type` is the type of plot to generate. One of `heatmap`,
    `box`, `violin` or `merged_heatmap` (all subsets in one plot).<br />
    * `devpars` is a dictionary with the device parameters for the plot.<br />
    * Other arguments for `plotthis::Heatmap()`, `plotthis::BoxPlot()`
    or `plotthis::ViolinPlot()`, depending on the `plot_type`.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    Multiple cases for the analysis.<br />
    If you only have one case, you can specify the parameters directly to
    `envs.ntimes`, `envs.subset_by`, `envs.group_by`, `envs.group1`,
    `envs.group2`, and `envs.plots`. The name of the case will be
    `envs.subset_by`.<br />
    If you have multiple cases, you can specify the parameters for each case
    in a dictionary. The keys will be the names of the cases and the values
    will be dictionaries with the parameters for each case, where the values
    will be inherited from `envs.ntimes`, `envs.subset_by`, `envs.group_by`,
    `envs.group1`, `envs.group2`, and `envs.plots`.<br />

