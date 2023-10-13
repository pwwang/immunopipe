# ScrnaMetabolicLandscape

Metabolic landscape analysis for scRNA-seq data

An abstract from <https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape>.

See also <https://pwwang.github.io/biopipen/pipelines/scrna_metabolic/>.

This is a [group of processes][1] to analyze the metabolic landscape of single cell RNA-seq data.

It collects a set of processes and owns a set of arguments. These arguments could either preset the default values for the processes or define the relationships between the processes.

## Processes

The processes in this group implement part of the pipeline below from the original paper.

![scrna_metabolic_pipeline](https://raw.githubusercontent.com/LocasaleLab/Single-Cell-Metabolic-Landscape/master/pipeline.png)

The data preparation, preprocessing and clustering already done by other processes of this pipeline. The processes in this group are used to analyze the metabolic landscape of the data.

- [`MetabolicInput`](./MetabolicInput.md): Input for the metabolic pathway analysis pipeline for scRNA-seq data
- [`MetabolicExprImpution`](./MetabolicExprImpution.md): Impute the missing values in the  expression data
- [`MetabolicPathwayActivity`](./MetabolicPathwayActivity.md): Calculate the pathway activities for each group
- [`MetabolicPathwayHeterogeneity`](./MetabolicPathwayHeterogeneity.md): Calculate the pathway heterogeneity
- [`MetabolicFeatures`](./MetabolicFeatures.md): Inter-subset metabolic features - Enrichment analysis in details
- [`MetabolicFeaturesIntraSubset`](./MetabolicFeaturesIntraSubset.md): Intra-subset metabolic features

## Group arguments

- `noimpute` (`flag`): Whether to do imputation for the dropouts.
    If `False`, the values will be left as is.
- `gmtfile`: The GMT file with the metabolic pathways. The gene names should
    match the gene names in the gene list in `RNAData` or the `Seurat` object
- `grouping`: It defines the basic groups to investigate the metabolic activity
    Typically the clusters.
- `grouping_prefix`: Working as a prefix to group names
    For example, if we have `grouping_prefix = "cluster"` and
    we have `1` and `2` in the `grouping` column, the groups
    will be named as `cluster_1` and `cluster_2`
- `subsetting` (`type=auto`): How do we subset the data. Other columns in the
    metadata to do comparisons. For example, `"TimePoint"` or
    `["TimePoint", "Response"]`
- `subsetting_prefix` (type=auto): Working as a prefix to subset names
    For example, if we have `subsetting_prefix = "timepoint"` and
    we have `pre` and `post` in the `subsetting` column, the subsets
    will be named as `timepoint_pre` and `timepoint_post`
    If `subsetting` is a list, this should also be a same-length
    list. If a single string is given, it will be repeated to a list
    with the same length as `subsetting`
- `subsetting_comparison` (type=json): What kind of comparisons are we
    doing to compare cells from different subsets.
    It should be dict with keys as the names of the comparisons and
    values as the 2 comparison groups from the `subsetting` column.
    For example, if we have `pre` and `post` in the `subsetting` column,
    we could have
    `subsetting_comparison = {"pre_vs_post": ["post", "pre"]}`
    The second group will be the control group in the comparison.
    If we also have `1`, `2` and `3` in the `grouping` column,
    by default, the comparisons are done within each subset for
    each group. For example, for group `1`, groups `2` and `3`
    will be used as control, and for group `2`, groups `1` and `3`
    will be used as control, and for group `3`, groups `1` and `2`
    will be used as control. It is similar to `Seurat::FindMarkers`
    procedure. With this option, the comparisons are also done to
    compare cells from different subsets within each group. With the
    example above, we will have `pre_vs_post` comparisons within
    each group.
    If `subsetting` is a list, this must be a list of dicts with the
    same length.
- `mutaters` (type=json): Add new columns to the metadata for
    grouping/subsetting.
    They are passed to `sobj@meta.data |> mutate(...)`. For example,
    `{"timepoint": "if_else(treatment == 'control', 'pre', 'post')"}`
    will add a new column `timepoint` to the metadata with values of
    `pre` and `post` based on the `treatment` column.
- `ncores` (type=int): Number of cores to use for parallelization for
    each process

## Reference

- [Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.][2]

[1]: https://pwwang.github.io/pipen/proc-group/
[2]: https://www.nature.com/articles/s41467-019-11738-0
