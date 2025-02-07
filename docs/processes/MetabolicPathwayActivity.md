# MetabolicPathwayActivity

This process calculates the pathway activities in different groups and subsets.

The cells are first grouped by subsets and then the metabolic activities are
examined for each groups in different subsets.<br />

For each subset, a heatmap and a violin plot will be generated.<br />
The heatmap shows the pathway activities for each group and each metabolic pathway

![MetabolicPathwayActivity_heatmap](../latest/processes/images/MetabolicPathwayActivity_heatmap.png){: width="80%"}

The violin plot shows the distribution of the pathway activities for each group

![MetabolicPathwayActivity_violin](../latest/processes/images/MetabolicPathwayActivity_violin.png){: width="45%"}

## Input

- `sobjfile`:

## Output

- `outdir`: *Default: `{{in.sobjfile | stem}}.pathwayactivity`*. <br />

## Environment Variables

- `ntimes` *(`type=int`)*: *Default: `5000`*. <br />
    Number of times to do the permutation
- `ncores` *(`type=int;pgarg`)*: *Default: `1`*. <br />
    Number of cores to use for parallelization
    Defaults to `ScrnaMetabolicLandscape.ncores`
- `heatmap_devpars` *(`ns`)*:
    Device parameters for the heatmap
    - `width` *(`type=int`)*:
        Width of the heatmap
    - `height` *(`type=int`)*:
        Height of the heatmap
    - `res` *(`type=int`)*:
        Resolution of the heatmap
- `violin_devpars` *(`ns`)*:
    Device parameters for the violin plot
    - `width` *(`type=int`)*:
        Width of the violin plot
    - `height` *(`type=int`)*:
        Height of the violin plot
    - `res` *(`type=int`)*:
        Resolution of the violin plot
- `gmtfile` *(`pgarg`)*:
    The GMT file with the metabolic pathways.<br />
    Defaults to `ScrnaMetabolicLandscape.gmtfile`
- `grouping` *(`type=auto;pgarg;readonly`)*:
    Defines the basic groups to
    investigate the metabolic activity, typically the clusters.<br />
    Defaults to `ScrnaMetabolicLandscape.grouping`
- `grouping_prefix` *(`type=auto;pgarg;readonly`)*: *Default: `""`*. <br />
    Working as a prefix to group
    names. For example, if we have `grouping_prefix = "cluster"` and
    we have `1` and `2` in the `grouping` column, the groups
    will be named as `cluster_1` and `cluster_2`.<br />
    Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
- `subsetting` *(`type=auto;pgarg;readonly`)*:
    How do we subset the data. Other
    columns in the metadata to do comparisons. For example,
    `"TimePoint"` or `["TimePoint", "Response"]`.<br />
    Defaults to `ScrnaMetabolicLandscape.subsetting`
- `subsetting_prefix` *(`type=auto;pgarg;readonly`)*:
    Working as a prefix to
    subset names.<br />
    For example, if we have `subsetting_prefix = "timepoint"` and
    we have `pre` and `post` in the `subsetting` column, the subsets
    will be named as `timepoint_pre` and `timepoint_post`.<br />
    If `subsetting` is a list, then this should also be a
    same-length list. If a single string is given, it will be
    repeated to a list with the same length as `subsetting`.<br />
    Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

