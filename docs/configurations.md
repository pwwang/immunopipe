
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

### `MARKERS_FINDER`

- See: [https://pwwang.github.io/immunopipe/markers-finder/](https://pwwang.github.io/immunopipe/markers-finder/)

### `GENE_EXPR_INVESTIGATION_CLUSTERS`

- See: [https://pwwang.github.io/immunopipe/gene-expr-investigation-for-each-cluster/](https://pwwang.github.io/immunopipe/gene-expr-investigation-for-each-cluster/)

### `RADAR_PLOTS`

- See: [https://pwwang.github.io/immunopipe/radar/](https://pwwang.github.io/immunopipe/radar/)]

### `METABOLIC`

- See: [https://pwwang.github.io/immunopipe/metabolic/](https://pwwang.github.io/immunopipe/metabolic/)]


## An examplar configuration file

```toml
 # pipeline options
forks = 4

# plugin options
[plugin_opts]
report_forks = 4
args_hide = true

# process settings
[SampleInfo.in]
infile = ["sample.txt"]

[CloneResidency.envs]
subject = ["patient"]
group = "region"
order = ["Normal adjacent tissue", "Tumor"]

[ImmunarchAdvanced.envs]
div_methods = ["gini", "div"]

[MarkersFinder.envs]
ncores = 4

[MarkersForClustersOfAllCells.envs]
ncores = 4

[MarkersForClustersOfTCells.envs]
ncores = 4

[SeuratPreparing.envs]
ncores = 12

[Immunarch.envs.tracking_target]
L1_top_5 = {TOP = 5}
L2_top_5 = {TOP = 5}
L3_top_5 = {TOP = 5}
L4_top_5 = {TOP = 5}
L5_top_5 = {TOP = 5}
L6_top_5 = {TOP = 5}

[Immunarch.envs.tracking_samples]
L1_top_5 = ["LN1", "LT1"]
L2_top_5 = ["LN2", "LT2"]
L3_top_5 = ["LN3", "LT3"]
L4_top_5 = ["LN4", "LT4"]
L5_top_5 = ["LN5", "LT5"]
L6_top_5 = ["LN6", "LT6"]

[DimPlots.envs.cases.Ident_UMAP]
"group.by" = "ident"
reduction = "umap"

[CloneHeterogeneity.envs.cases.By_Quantile]
cut = "quantile"

[CloneHeterogeneity.envs.cases.By_Quantile.subsetting]
Tumor = 'region == "Tumor"'
Normal = 'region == "Normal adjacent tissue"'

[CloneHeterogeneity.envs.cases.By_Quantile.design]
Tumor_vs_Normal = ["Tumor", "Normal"]

# Module based settings

[[MARKERS_FINDER]]  # 0
name = "Markers_of_top_20_TCR_clones"
[MARKERS_FINDER.design]
LT1_vs_LN1 = ["LT1", "LN1"]
LT2_vs_LN2 = ["LT2", "LN2"]
LT3_vs_LN3 = ["LT3", "LN3"]
LT4_vs_LN4 = ["LT4", "LN4"]
LT5_vs_LN5 = ["LT5", "LN5"]
LT6_vs_LN6 = ["LT6", "LN6"]
[MARKERS_FINDER.meta]
MetaMarkers_Tumor = ["LT1", "LT2", "LT3", "LT4", "LT5", "LT6"]
MetaMarkers_Normal = ["LN1", "LN2", "LN3", "LN4", "LN5", "LN6"]
[MARKERS_FINDER.overlap]
SharedMarkers = ["LT1_vs_LN1", "LT2_vs_LN2", "LT3_vs_LN3", "LT4_vs_LN4", "LT5_vs_LN5", "LT6_vs_LN6"]
[MARKERS_FINDER.filters.LT1]
"by.meta" = { Sample = "LT1" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LT2]
"by.meta" = { Sample = "LT2" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LT3]
"by.meta" = { Sample = "LT3" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LT4]
"by.meta" = { Sample = "LT4" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LT5]
"by.meta" = { Sample = "LT5" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LT6]
"by.meta" = { Sample = "LT6" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LN1]
"by.meta" = { Sample = "LN1" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LN2]
"by.meta" = { Sample = "LN2" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LN3]
"by.meta" = { Sample = "LN3" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LN4]
"by.meta" = { Sample = "LN4" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LN5]
"by.meta" = { Sample = "LN5" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}
[MARKERS_FINDER.filters.LN6]
"by.meta" = { Sample = "LN6" }
"by.count" = { ORDER = 1, filter = "CDR3.aa %in% CDR3.aa[1:10]"}

# Gene expression investigation for clusters
[[GENE_EXPR_INVESTIGATION_CLUSTERS]] # 0
genefile = "genes.txt"
subset = "region == 'Tumor'"
name = "[Tumor][boxplots] Gene expressions of T-Cell clusters"
groupby = "seurat_clusters"
plots = { boxplot = {ncol = 3, res = 100, width = 1200, height = 1200} }

[[GENE_EXPR_INVESTIGATION_CLUSTERS]] # 1
genefile = "hmgenes.txt"
subset = "region == 'Tumor'"
name = "[Tumor][heatmap] Gene expressions of T-Cell clusters"
groupby = "seurat_clusters"
plots = { heatmap = {res = 100, width = 1200, height = 1200} }

[[GENE_EXPR_INVESTIGATION_CLUSTERS]] # 0
genefile = "genes.txt"
subset = "region == 'Normal adjacent tissue'"
name = "[Normal][boxplots] Gene expressions of T-Cell clusters"
groupby = "seurat_clusters"
plots = { boxplot = {ncol = 3, res = 100, width = 1200, height = 1200} }

[[GENE_EXPR_INVESTIGATION_CLUSTERS]] # 1
genefile = "hmgenes.txt"
subset = "region == 'Normal adjacent tissue'"
name = "[Normal][heatmap] Gene expressions of T-Cell clusters"
groupby = "seurat_clusters"
plots = { heatmap = {res = 100, width = 1200, height = 1200} }

[[RADAR_PLOTS]]
name = "Cell proportion"
breaks = [0, 15, 30]
[RADAR_PLOTS.filters.BM]
"by.meta" = { region = "Tumor" }
[RADAR_PLOTS.filters.PB]
"by.meta" = { region = "Normal adjacent tissue"}

[METABOLIC]
gmtfile = "KEGG_metabolism.gmt"

[[METABOLIC.cases]]

[METABOLIC.cases.grouping]
groupby = "idents"

[METABOLIC.cases.subsetting]
groupby = "region"
alias = "SampleType"

[METABOLIC.cases.design]
Tumor_vs_Normal = ["Tumor", "Normal adjacent tissue"]

[MetabolicPathwayActivity.envs.heatmap_devpars]
width = 1200
height = 1000

```