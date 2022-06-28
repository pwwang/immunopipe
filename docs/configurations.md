
## Arguments from command line

`immunopipe` is written using `pipen`. The pipeline is composed of a number processes. Each process' input (start processes), envs and other configurations can be passed from command line.

For example, to run V-J usage analysis on 8 cores:
```shell
immunopipe --VJUsage.envs.ncores 8 ...
```

## Arguments from configuration file

As the argument parsing for pipen is supported by `pipen-args`, which allows us to pass arguments from both command line and a configuration file.

To specify the configuration file:
```shell
immunopipe --config config.toml ...
```

The configuration file is in [TOML](https://github.com/toml-lang/toml) format.

You can specify the same argument from both command line and the configuration file, but the one from the command line has higher priority.

To specify the arguments, for example, `--VJUsage.envs.cores=8` from the configuration file:

```toml
[VJUsage.envs]
cores = 8
```

You can also specify some arguments for the entire pipeline, for example, the scheduler of the pipeline:

```toml
scheduler = "sge"

[VJUsage.envs]
cores = 8
```

## Exploration of all configurable items for the pipeline or for a specific process

To see all available configuration items for the pipeline, run:

```shell
immunopipe --full
```

For a process, for example, for `VJUsage`:

```shell
pipen run tcr VJUsage --full
```

Or you can check the docstring of the process in `biopipen`'s source code.

## Configurations for other modules

### `METABOLIC`

- See: [https://pwwang.github.io/immunopipe/metabolic/](https://pwwang.github.io/immunopipe/metabolic/)]


## An examplar configuration file

```toml
# Full configurations for immunopipe
# ==================================
# # You can also run `immunopipe --config config.toml --full`
# # to see full configuration options from command line.

# Pipeline options
# ----------------
# # See: https://pwwang.github.io/pipen/configurations
# profile = 'default'
# outdir = './output'
# loglevel = 'info'
# cache = true
# dirsig = 1
# error_strategy = 'ignore'
# num_retries = 3
forks = 4
# submission_batch = 8
# name = 'immunopipe'
# scheduler = 'local'
# scheduler_opts = {}
# plugins = []
# template_opts = {}
[plugin_opts]
report_forks = 4
args_hide = true

# General process options
# -----------------------
# # See: https://pwwang.github.io/pipen/defining-proc/#process-configurations-and-proc-class-variables
# [<process_name>]
# forks = 1
# # ... other general process options

# Report options
# --------------
# # See: https://pwwang.github.io/pipen-report/configurations
# [<process_name>.plugin_opts]
# report_toc = false
# # ... other report options

# Process input data
# ------------------
# # Only start process can have input data specified
# # The only start process of immunopipe is 'SampleInfo'
# # See: https://pwwang.github.io/immunopipe/input
[SampleInfo.in]
infile = ["sample.txt"]

# Process environments
# --------------------
# # [<process_name>.envs]
# # env_name = env_value

[Immunarch.envs]
# # Groupings to show clonotype volume (sizes) Default: {}
# # Multiple groups supported, for example: volume_by = {0: "Status", 1: ["Status", "Sex"]}
# # Or label the groups: `volume_by = {"Status": "Status", "Status_Sex": ["Status", "Sex"]}`
# # If a list or a single variable is given, it will be changed into {"Status": "Status"}
# volume_by = "Status"  # or multiple groups: volume_by = {"Status_Sex": ["Status", "Sex"]}
# # Groupings to show CDR3 length of both aa and nt Default: {}
# len_by = "Status"
# # Groupings to show clonotype counts per sample
# count_by = "Status"
# # `.head` arguments of Immunarch::repClonoality()
# top_clone_marks = [10, 100, 1000, 3000, 10000]
# # Groupings when visualize top clones
# top_clone_by = "Status"
# # `.bound` arguments of `repClonoality()`
# rare_clone_marks = [1, 3, 10, 30, 100]
# # Groupings when visualize rare clones
# rare_clone_by = {}
# # `.clone.types` arguments of `repClonoality()`
# hom_clone_marks = {Small = 0.0001, Medium = 0.001, Large = 0.01, Hyperexpanded = 1}
# # Groupings when visualize homeo clones
# hom_clone_by = {}
# # The methods used for `repOverlap()`, each will generate a heatmap.
# overlap_methods = ["public"]
# # Plot the samples with these dimension reduction methods
# overlap_redim = ['tsne', 'mds']
# # Groupings to show gene usages
# gu_by = {}
# # How many top (ranked by total usage across samples) genes to show in the plots
# gu_top = 30
# # Controls how the data is going to be preprocessed and analysed.
# # Some of js, cor, cosine, pca, mds, and tsne
# gua_methods = ['js', 'cor']
# # .quant and .col for spectratype() for each sample
# spect = [{quant = 'id', col = 'nt'}, {quant = count, col = 'aa+v'}]
# # Methods to calculate diversities
# div_methods = ['div', 'gini.simp']
# # Groupings to show sample diversities
# div_by = {}
# # Groupings to show rarefactions
# raref_by = {}
# # The target and samples to track. You can do multiple trackings.
# # To do that, you need to specify a key for each tracking.
# # It will use the target and samples under the same key.
# # If samples from `tracking_samples` cannot be found, all samples will be used
# # Other than the target supported by immunarch, you can also specify top shared clones. For example:
[Immunarch.envs.tracking_target]
L1_top_5 = { TOP = 5 }
L2_top_5 = { TOP = 5 }
L3_top_5 = { TOP = 5 }
L4_top_5 = { TOP = 5 }
L5_top_5 = { TOP = 5 }
L6_top_5 = { TOP = 5 }
[Immunarch.envs.tracking_samples]
L1_top_5 = ["LN1", "LT1"]
L2_top_5 = ["LN2", "LT2"]
L3_top_5 = ["LN3", "LT3"]
L4_top_5 = ["LN4", "LT4"]
L5_top_5 = ["LN5", "LT5"]
L6_top_5 = ["LN6", "LT6"]
# # Arguments for kmer analysis.
# # Keys are the K of mers. Values are parameters:
# # - `head` specifies # of the most abundant kmers to visualise.
# # - `position`: positions of bars: `stack`, `dodge` and `fill`
# # - `log`: log-transformation of y-axis
# # - `motif`: Method for motif analysis
# # There can be multiple `head`s and `motif`s.
# # If you do want multiple parameter sets for the same K, You can use
# # a float number as the K. For example: `5.1` for K `5`.
# kmers = {5 = {head = 10, position = "stack", log = false, motif = "self"}}

[SeuratPreparing.envs]
# # Number of cores to use
ncores = 4
# # Seurat routine: sct or integrate
# # sct routine - SelectIntegrationFeatures -> PrepSCTIntegration -> FindIntegrationAnchors -> IntegrateData -> RunPCA -> RunUMAP
# # integrate routine - NormalizeData -> FindVariableFeatures -> FindIntegrationAnchors -> ScaleData -> RunPCA -> RunUMAP
# routine = "integrate"
# # Arguments for Seurat functions in the routine. For example:
IntegrateData = {"k.weight" = 20}

# [SeuratClusteringOfAllCells.envs]
# # Arguments for Seurat::FindNeighbors
# FindNeighbors = {}
# # Arguments for Seurat::FindClusters
# FindClusters = {resollution: 0.8}

[MarkersForClustersOfAllCells.envs]
# # Number of cores to use
ncores = 4
# # The dbs to do enrichment analysis for significant markers
# # See: https://maayanlab.cloud/Enrichr/#libraries
dbs = ["KEGG_2021_Human"]

[SelectTCells.envs]
# Clonotype percentage threshold to filter cells
tcell_filter = "Clonotype_pct > 0.1"
# Gene expression used to filter cells
indicator_gene = "CD3E"

# [SeuratClusteringOfTCells.envs]
# # Same as SeuratClusteringOfAllCells.envs

[MarkersForClustersOfTCells.envs]
# # Same as MarkersForClustersOfAllCells
ncores = 4
dbs = ["KEGG_2021_Human"]

[SeuratClusterStats.envs]
# # Statistics about Seurat clusters
[SeuratClusterStats.envs.stats]
# # The statistics to plot
# # - nCells - Number of cells for each cluster
# # - nCellsPerSample - Number of cells per sample for each cluster
# # - percCellsPerSample - Percentage of cells per sample for each cluster
nCells = {res = 100, height = 1000, width = 1000}
nCellsPerSample = {res = 100, height = 1000, width = 1000}
percCellsPerSample = {res = 100, height = 1000, width = 1000}
# # The expression values to plot
[SeuratClusterStats.envs.exprs]
# # ridgeplots - The ridge plots for the gene expressions.
# # See `?Seurat::RidgePlot`.
# # vlnplots - Violin plots for the gene expressions.
# # See `?Seurat::VlnPlot`. You can have `boxplot` key to add
# # `geom_boxplot()` to the violin plots
# # featureplots - The feature plots for the gene expressions.
# # See `?Seurat::FeaturePlot`.
# # dotplot - Dot plots for the gene expressions.
# # See `?Seurat::DotPlot`.
# # heatmap - Heatmap for the gene expressions.
# # See `?Seurat::DoHeatmap`. You can specify `average=True` to plot on
# # the average of the expressions.
# # All the above with `devpars` to define the output figures
# # and `plus` to add elements to the `ggplot` object.
# # You can have `subset` to subset the data. Multiple cases can be
# # distinguished by `ridgeplots` and `ridgeplots.1`
"ridgeplots.1"  = {title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"ridgeplots.2"  = {title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"vlnplots.1"  = {boxplot = {}, "pt.size" = 0, title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"vlnplots.2"  = {boxplot = {}, "pt.size" = 0, title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"featureplots.1"  = {title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"featureplots.2"  = {title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"dotplot.1"  = {plus = "RotatedAxis()", title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"dotplot.2"  = {plus = "RotatedAxis()", title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
"heatmap.1"  = {downsample = "average", title = "Gene expressions in Tumor", subset = "region == 'Tumor'"}
"heatmap.2"  = {downsample = "average", title = "Gene expressions in Normal", subset = "region == 'Normal adjacent tissue'"}
# # The dimensional reduction plots
[SeuratClusterStats.envs.dimplots]
# # `<case>` - The case to plot. Keys are the arguments for `Seurat::Dimplot()`, add `devpars`.
Ident = {"group.by" = "ident", devpars = {res = 100, width = 1000, height = 1000}}

[TCRClustering.envs]
# # The tool used to do the clustering, either GIANA or ClusTCR
tool = "ClusTCR"
# # Arguments for GIANA/ClusTCR
# # For GIANA, they will be passed to `python GIAna.py <args>`
# # For ClusTCR, they will be passed to `clustcr.Clustering(<args>)`
# args = {}

[TCRClusteringStats.envs]
# # numbers_on_heatmap: Whether to show the numbers on the heatmap
# # heatmap_meta: The metadata to show on the heatmap
# # grouping: The groups to investigate the shared clusters
shared_clusters = { numbers_on_heatmap = true, heatmap_meta = [], grouping = "null" }
# # Sample diversity using TCR clusters instead of clones
# # keys are the methods and values, currently, `by` to plot
# # the diversities by groups
sample_diversity = { gini = { by = ["region"] } }

[CellsDistribution.envs]
# # Show cells/clone/TCR cluster distribution in Seurat clusters
# # Cases (subsets/specific clones/clusters/etc)
[CellsDistribution.envs.cases.Tumor_vs_Normal]
# # Columns of the pie charts
[CellsDistribution.envs.cases.Tumor_vs_Normal.group]
by = "Region"
order = ["Tumor", "Normal"]
# # You can also create new columns in the data
[CellsDistribution.envs.cases.Tumor_vs_Normal.group.mutaters]
Region = "if_else (region == 'Tumor', 'Tumor', 'Normal')"
# # Rows of the pie charts
[CellsDistribution.envs.cases.Tumor_vs_Normal.cells]
# # Select the cells/clone/clusters by?
by = "TCR_Cluster"
# # How many?
n = 10
# # order by? Will be passed to `dplyr::arrange()`
# # Available size columns to use: .CloneSize, .CloneGroupSize, .CloneGroupClusterSize
orderby = "desc(.CloneSize)"

[CloneResidency.envs]
# # How are TCR clones distributed in different groups
# # The key of subject in metadata. The clone residency will be examined for this subject/patient
subject = ["patient"]
# # The key of group in metadata. This usually marks the samples that you want to compare. For example, Tumor vs Normal,
# # post-treatment vs baseline. It doesn't have to be 2 groups always. If there are more than 3
# # groups, instead of venn diagram, upset plots will be used.
group = "region"
# # The order of the values in `group`. Early-ordered group will be used as x-axis in scatter plots
# # If there are more than 2 groups, for example, [A, B, C], the scatter plots will be drawn for pairs: B ~ A, C ~ B.
order = ["Normal adjacent tissue", "Tumor"]

[CloneHeterogeneity.envs]
# # Clone investigation in Seurat clusters
# # Cases (subsets/specific clones/clusters/etc)
[CloneHeterogeneity.envs.cases.By_Quantile]
# # How to cut the clones by sizes
cut = "quantile"
# # Subsets to use
[CloneHeterogeneity.envs.cases.By_Quantile.subsetting]
Tumor = 'region == "Tumor"'
Normal = 'region == "Normal adjacent tissue"'
# # Designed comparisons
[CloneHeterogeneity.envs.cases.By_Quantile.design]
Tumor_vs_Normal = ["Tumor", "Normal"]

[RadarPlots.envs]
# # Radar plots on stats for each Seurat cluster
# # Cases (subsets/specific clones/clusters/etc)
[RadarPlots.envs.cases.Tumor_vs_Normal]
# # Add new columns to the meta.data
mutaters = {Group = "if_else (region == 'Tumor', 'Tumor', 'Normal')"}
# # Which column to use to separate the cells in different groups
by = "Group"
# # The order of the values in `by`.
order = ["Tumor", "Normal"]
# # breaks of the radar plots
breaks = [0, 15, 30]
# # Direction to calculate the percentages
# # inter-cluster: the percentage of the cells in all groups
# # intra-cluster: the percentage of the cells in all clusters
direction = "intra-cluster"
# # Whether the percentages are calculated before or after filtering out the NAs
perc_with_na = false

[ScFGSEA.envs]
# # Do GSEA on two groups of cells
# # Like seudo-bulk GSEA
# # gmtfile with pathways or gene sets
gmtfile = "MSigDB_Hallmark_v7.5.1.gmt"
# # One could also use placeholders for the cases.
# # Currently only cluster is supported. One could use `{cluster}` or
# # `{ident}` to denote the clusters.
[ScFGSEA.envs.cases.Tumor_vs_Normal]
# # Add columns to meta data
mutaters = {Group = "if_else (region == 'Tumor', 'Tumor', 'Normal')"}
# # The first group of cells
"ident.1" = "Tumor"
# # The second group of cells
"ident.2" = "Normal"
# # From which column?
"group.by" = "Group"

[MarkersFinderClones.envs]
# # Do marker finding for clones
# # Similar as `ScFGSEA.envs.cases`. "ident.2" is optional. In such a case,
# # The rest of the cells are used as the second group.
[MarkersFinderClones.envs.cases.Top20_Clones_Tumor]
"ident.1" = "Top20_Clones"
"ident.2" = "Rest_Clones"
"group.by" = "Top20_Clones"
[MarkersFinderClones.envs.cases.Top20_Clones_Tumor.mutaters]
CloneSize = "group_by(meta, CDR3.aa) |> mutate(n=if_else(region != 'Tumor', NA_integer_, n())) |> pull('n')"
Top20_Clones = """if_else (
  CloneSize %in% sort(CloneSize, decreasing=TRUE)[1:20],
  'Top20_Clones',
  'Rest_Clones'
)"""
[MarkersFinderClones.envs.cases.Top20_Clones_Normal]
"ident.1" = "Top20_Clones"
"ident.2" = "Rest_Clones"
"group.by" = "Top20_Clones"
[MarkersFinderClones.envs.cases.Top20_Clones_Normal.mutaters]
CloneSize = "group_by(meta, CDR3.aa) |> mutate(n=if_else(region != 'Normal adjacent tissue', NA_integer_, n())) |> pull('n')"
Top20_Clones = """if_else (
  CloneSize %in% sort(CloneSize, decreasing=TRUE)[1:20],
  'Top20_Clones',
  'Rest_Clones'
)"""

# Module options
# --------------

[METABOLIC]
# # The metabolic pathway in gmt file
gmtfile = "KEGG_metabolism.gmt"
[[METABOLIC.cases]]
# # The name of the case
name = "Case1"
# # How do we group the cells?
[METABOLIC.cases.grouping]
# # Add new columns to the meta.data
# mutaters = {}
# # The columns to group the cells
groupby = "seurat_clusters"
# # How do we subset the data. The imputation will be done in each subset separately
[METABOLIC.cases.subsetting]
# # Add new columns to the meta.data
# mutaters = {}
# # The columns to subset the data
groupby = "region"
# # The alias of the subset working as a prefix to subset names
alias = "SampleType"
# # What kind of comparisons are we doing?
# # It should be the values of subsetting `groupby`s
[METABOLIC.cases.design]
Tumor_vs_Normal = ["Tumor", "Normal adjacent tissue"]

# # More cases
# [[METABOLIC.cases]]
# name = "Case2"
# ...

# Options for module processes
# ----------------------------

[MetabolicPathwayActivity.envs.heatmap_devpars]
width = 1200
height = 1000

```