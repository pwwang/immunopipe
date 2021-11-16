
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

### `DIM_PLOTS`

- See: [https://pwwang.github.io/immunopipe/dimplots/](https://pwwang.github.io/immunopipe/dimplots/)]

### `RADAR_PLOTS`

- See: [https://pwwang.github.io/immunopipe/radar/](https://pwwang.github.io/immunopipe/radar/)]


## An examplar configuration file

```toml
# pipeline options
forks = 4

# plugin options
[plugin_opts]
report_forks = 4

# process settings
[SampleInfo.in]
infile = ["samples.txt"]

[GeneList.in]
metafile = ["genes.txt", "hmgenes.txt"]

[CloneResidency.envs]
subject = ["Patient"]
group = "Source"
order = ["PB", "BM"]

[ImmunarchAdvanced.envs]
div_methods = ["gini", "div"]

[MarkersFinder.envs]
ncores = 4

[MarkersForClustersOfAllCells.envs]
ncores = 32

[MarkersForClustersOfTCells.envs]
ncores = 32

[SeuratPreparing.envs]
ncores = 12

[ImmunarchAdvanced.envs.tracking_target]
MM005_top_10 = {TOP = 10}
MM006_top_10 = {TOP = 10}

[ImmunarchAdvanced.envs.tracking_samples]
MM005_top_10 = ["MM005BM-earlier", "MM005BM-postr", "MM005WBC-earlier", "MM005WBC-postr"]
MM006_top_10 = ["MM006BM-pre", "MM006BM-postr", "MM006WBC-pre", "MM006WBC-postr"]

# Module based settings

## Pre-CART BM
[[MARKERS_FINDER.filters]]  # 0
name = "[pre-CART BM DR] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "BM", TimePoint = "Pre-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "BM", TimePoint = "Pre-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 1
name = "[pre-CART BM PD] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "BM", TimePoint = "Pre-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "BM", TimePoint = "Pre-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 2
name = "[BM CNTRL] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "BM", Response = "CNTRL" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "BM", Response = "CNTRL" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

## Pre-CART PB
[[MARKERS_FINDER.filters]]  # 3
name = "[pre-CART PB DR] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "PB", TimePoint = "Pre-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "PB", TimePoint = "Pre-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 4
name = "[pre-CART PB PD] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "PB", TimePoint = "Pre-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "PB", TimePoint = "Pre-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 5
name = "[PB CNTRL] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "PB", Response = "CNTRL" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "PB", Response = "CNTRL" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

## Post-CART BM
[[MARKERS_FINDER.filters]]  # 6
name = "[Post-CART BM DR] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "BM", TimePoint = "Post-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "BM", TimePoint = "Post-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 7
name = "[Post-CART BM PD] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "BM", TimePoint = "Post-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "BM", TimePoint = "Post-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

## Post-CART PB
[[MARKERS_FINDER.filters]]  # 8
name = "[Post-CART PB DR] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "PB", TimePoint = "Post-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "PB", TimePoint = "Post-CART",  Response = "DR" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 9
name = "[Post-CART PB PD] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "PB", TimePoint = "Post-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "PB", TimePoint = "Post-CART",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

## Early-BMT BM
# [[MARKERS_FINDER.filters]]  # 10
# name = "[Early-BMT BM DR] Markers of top 20 TCR clones"
# [MARKERS_FINDER.filters."ident.1"]
# "by.meta" = { Source = "BM", TimePoint = "Early-BMT",  Response = "DR" }
# "by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
# [MARKERS_FINDER.filters."ident.2"]
# "by.meta" = { Source = "BM", TimePoint = "Early-BMT",  Response = "DR" }
# "by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 11
name = "[Early-BMT BM PD] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "BM", TimePoint = "Early-BMT",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "BM", TimePoint = "Early-BMT",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

## Early-BMT PB
# [[MARKERS_FINDER.filters]]  # 12
# name = "[Early-BMT PB DR] Markers of top 20 TCR clones"
# [MARKERS_FINDER.filters."ident.1"]
# "by.meta" = { Source = "PB", TimePoint = "Early-BMT",  Response = "DR" }
# "by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
# [MARKERS_FINDER.filters."ident.2"]
# "by.meta" = { Source = "PB", TimePoint = "Early-BMT",  Response = "DR" }
# "by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

[[MARKERS_FINDER.filters]]  # 13
name = "[Early-BMT PB PD] Markers of top 20 TCR clones"
[MARKERS_FINDER.filters."ident.1"]
"by.meta" = { Source = "PB", TimePoint = "Early-BMT",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "TOTAL %in% TOTAL[1:20]"}
[MARKERS_FINDER.filters."ident.2"]
"by.meta" = { Source = "PB", TimePoint = "Early-BMT",  Response = "PD" }
"by.count" = { ORDER = 1,  filter = "!TOTAL %in% TOTAL[1:20]" }

# Gene expression investigation for clusters
[[GENE_EXPR_INVESTIGATION_CLUSTERS]] # 0
name = "[Control samples] Gene expressions of T-Cell clusters"
target = [
    "Cone051520A", "Cone051520B", "Cone051920", "Cone052720", "Cone061920",
    "WIG022221BM", "WIG022221WBC", "WIG030821BM", "WIG030821WBC", "WIG021721BM",
    "WIG021721WBC",
]
plots = { boxplot = {ncol = 3, res = 100, width = 1200, height = 3200} }

[[GENE_EXPR_INVESTIGATION_CLUSTERS]] # 1
name = "[Case samples] Gene expressions of T-Cell clusters"
target = [
    "MM001BM-earlier", "MM001BM-postr", "MM002BM-earlier", "MM003BM-earlier",
    "MM003WBC-earlier", "MM003BM-pre", "MM004BM-pre", "MM005BM-earlier",
    "MM005WBC-earlier", "MM005BM-postr", "MM005WBC-postr", "MM006BM-pre",
    "MM006WBC-pre", "MM006BM-postr", "MM006WBC-postr", "MM007BM-pre",
    "MM007BM-post", "MM007WBC-post", "MM008BM-pre", "MM009BM-pre",
    "MM009WBC-pre", "MM010BM-pre", "MM010WBC-pre", "MM011BM-earlier",
    "MM011WBC-earlier", "MM011BM-postr", "MM012BM-pre", "MM012WBC-pre",
    "MM013BM-pre", "MM013BM-post", "MM013WBC-post", "MM014BM-pre",
    "MM014WBC-pre", "MM015BM-post", "MM015WBC-post",
]
plots = { boxplot = {ncol = 3, res = 100, width = 1200, height = 3200} }

[[DIM_PLOTS]]
name = "Dimension redduction plot for clusters"
reduction = "umap"
"by.ident" = true

[[RADAR_PLOTS]]
name = "Cell proportion of control samples"
breaks = [0, 15, 30]
[RADAR_PLOTS.filters.BM]
"by.meta" = { Source = "BM", Status = "CNTRL" }
[RADAR_PLOTS.filters.PB]
"by.meta" = { Source = "PB", Status = "CNTRL" }

[[RADAR_PLOTS]]
name = "Cell proportion of pre-CART early relapse samples"
breaks = [0, 15, 30]
[RADAR_PLOTS.filters.BM]
"by.meta" = { Source = "BM", Status = "Pre", PFS_12Mo = "PFS_12Mo-Y" }
[RADAR_PLOTS.filters.PB]
"by.meta" = { Source = "PB", Status = "Pre", PFS_12Mo = "PFS_12Mo-Y" }

[[RADAR_PLOTS]]
name = "Cell proportion of pre-CART durable response samples"
breaks = [0, 15, 30]
[RADAR_PLOTS.filters.BM]
"by.meta" = { Source = "BM", Status = "Pre", PFS_12Mo = "PFS_12Mo-N" }
[RADAR_PLOTS.filters.PB]
"by.meta" = { Source = "PB", Status = "Pre", PFS_12Mo = "PFS_12Mo-N" }

[[RADAR_PLOTS]]
name = "Cell proportion of post-CART early relapse samples"
breaks = [0, 15, 30]
[RADAR_PLOTS.filters.BM]
"by.meta" = { Source = "BM", Status = "Post", PFS_12Mo = "PFS_12Mo-Y" }
[RADAR_PLOTS.filters.PB]
"by.meta" = { Source = "PB", Status = "Post", PFS_12Mo = "PFS_12Mo-Y" }

[[RADAR_PLOTS]]
name = "Cell proportion of post-CART durable response samples"
breaks = [0, 15, 30]
[RADAR_PLOTS.filters.BM]
"by.meta" = { Source = "BM", Status = "Post", PFS_12Mo = "PFS_12Mo-N" }
[RADAR_PLOTS.filters.PB]
"by.meta" = { Source = "PB", Status = "Post", PFS_12Mo = "PFS_12Mo-N" }

```