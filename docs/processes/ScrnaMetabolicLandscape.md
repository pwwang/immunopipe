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
- [`MetabolicExprImputation`](./MetabolicExprImputation.md): Impute the missing values in the  expression data
- [`MetabolicPathwayActivity`](./MetabolicPathwayActivity.md): Calculate the pathway activities for each group
- [`MetabolicPathwayHeterogeneity`](./MetabolicPathwayHeterogeneity.md): Calculate the pathway heterogeneity
- [`MetabolicFeatures`](./MetabolicFeatures.md): Inter-subset metabolic features - Enrichment analysis in details

## Group arguments

- `noimpute` (`flag`): Whether to do imputation for the dropouts.
    If `False`, the values will be left as is.
- `gmtfile`: The GMT file with the metabolic pathways. The gene names should
    match the gene names in the gene list in RNAData or
    the Seurat object.
    You can also provide a URL to the GMT file.
    For example, from
    <https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/>.
- `subset_by` (pgarg;readonly): Subset the data by the given column in the
    metadata. For example, `Response`.
    `NA` values will be removed in this column.
    If None, the data will not be subsetted.
- `group_by` (pgarg;readonly): Group the data by the given column in the
    metadata. For example, `cluster`.
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