# Configurations

In this section, we are discussing how to configure the pipeline itself and some common issues we may encounter or need to be aware of to configure individual processes. For the configurations of each process in details, please refer to the individual process pages.

As mentioned in `pipen`'s [docs](https://pwwang.github.io/pipen/configurations/), the configurations to control the pipeline can be `~/.pipen.toml` and/or `./.pipen.toml`. You can set the default values in those files. Additionally, you can also pass a configuration file from command line, as described in [`Running the pipeline`](./running.md) section. You can override the default values in the configuration files by passing arguments from command line.

The configuration file is in [`toml`](https://toml.io) format. A schematic example is shown below:

```toml
forks = 4
# Other pipeline configurations

# process configurations
[TCellSelection]
forks = 2  # override the default value

# envs of the process
# e.g
[TCellSelection.envs]
indicator_genes = ["CD3D", "CD3E", "CD3G"]

# other processes
# [ProcessName]
# ...
```

/// Tip
In the individual process pages, we will list the `envs` of the process. For example,

- `indicator_genes (list)`: The genes to be used to select T cells.

This means that the environment variable `indicator_genes` should be set as follows:

```toml
[TCellSelection.envs]
indicator_genes = ["CD3D", "CD3E", "CD3G"]
```
///


## Pipeline configurations

There are pipeline level configurations and process level configurations. The pipeline level configurations are used to control the pipeline itself. The process level configurations set here are the default values for all the processes. You can override the default values for each process in the process level configurations.

You can check all available configuration items and more details [here](https://pwwang.github.io/pipen/configurations/). Here we only list some of the most important ones. The rest ones are not recommended to change unless you know what you are doing.

### Pipeline level configurations

- `name`: The name of the pipeline (Default: `"Immunopipe"`)
  - It will change the working directory to `./.pipen/<name>`, where the pipeline information and intermediate files will be stored.
  - It will also change the default output directory to `./<name>-output`.
- `outdir`: The output directory (Default: `"./<name>-output"`)
    - See also [`Output directory and working directory`](#output-directory-and-working-directory).
- `loglevel`: The logging level for the logger (Default: `"info"`)
- `plugin_opts`: The options for the plugins.
    - Following `pipen` plugins are installed with `immunopipe`. You may check the links for more details.
    - [`pipen-board`][1]: Visualizing configuration and running of pipen pipelines on the web.
    - [`pipen-verbose`][2]: Adding verbosal information in logs for pipen.
    - [`pipen-runinfo`][3]: Generating running information for jobs in pipen pipelines.
    - [`pipen-filters`][4]: Adding a set of useful filters for pipen templates.
    - [`pipen-args`][5]: Command line argument parser for pipen
    - [`pipen-annotate`][6]: Using docstring to annotate pipen processes.
    - [`pipen-report`][7]: Generating reports for pipen pipelines.
    - [`pipen-log2file`][8]: Logging to files for pipen pipelines.
    - [`pipen-cli-run`][9]: Running pipen processes/process groups from command line.
- `scheduler`: The scheduler to use (Default: `"local"`)
- `scheduler_opts`: The options for the scheduler.
    - `immunopipe` is implemented using `pipen`, which is backended by [`xqute`][10]. Supported schedulers and options are listed [here](https://github.com/pwwang/xqute).
    - See also [`How to run the pipeline on a cluster?`](./faq.md#how-to-run-the-pipeline-on-a-cluster) for more details.

### Output and working directory

The output directory is the directory where the final results are stored. The working directory is the directory where the pipeline information and intermediate files are stored. By default, the output directory is `./<name>-output` and the working directory is `./.pipen/<name>`.

The output of processes with final results will be stored in the output directory in sub-directories named after the processes. For example, the output of [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) will be stored in `./<outdir>/SeuratClusteringOfAllCells`. This is also a special subdirectory named `REPORTS` that contains the reports of the pipeline. By visiting the `<outdir>/REPORTS/index.html` with a web browser, you can check the reports of the pipeline.

You can change the output directory by setting `outdir` or `name` in the configuration file. For example, if you want to change the output directory to `./output`, you can set the configurations as follows:

```toml
outdir = "./output"
```

If you change the pipeline name:

```toml
name = "my-pipeline"
```

Then the output directory will be changed to `./my-pipeline-output`.

/// Note
If both `outdir` and `name` are set, `outdir` will be used.
///

You can do the similar thing to change the working directory. However, you are **NOT** recommended to change the working directory, especially if you are using [`pipen-board`][1]. This is because that the plugin scans `./.pipen/<name>` to get the information for the previous run of the pipeline. If you change the working directory, the plugin will not be able to find the information for the previous run.

/// Tip
What if you want to change the working directory anyway? The recommended way is to create a symbolic link to the working directory. For example, if you want to change the working directory to `/path/to/the/real/working/directory`, you can do:

```bash
ln -s /path/to/the/real/working/directory ./.pipen
```
///

/// Tip
You can also then debug the pipeline by inspecting the real scripts in the working directory that run for the jobs of each process at `./.pipen/<name>/<process-name>/<job-index>/job.script`.

You can also find the other information for the jobs at `./.pipen/<name>/<process-name>/<job-index>/`, including the stdout (`job.stdout`) and stderr (`job.stderr`) of the jobs, the exit code of the jobs (`job.rc`), etc.
///

### Process level configurations

- `cache`: Should we detect whether the jobs are cached. If `true`, the jobs will be skipped if the output files exist and newer than the input files. (Default: `true`)
- `error_strategy`: The strategy to handle the errors.
    - `halt`: Any failure will just halt the entire pipeline (default)
    - `ignore`: Ignore the error and keep running (assuming the job runs successfully anyway)
    - `retry`: Retry to job running.
      After `num_retries` times of retrying, if the job is still failing, halt the pipeline.

- `num_retries`: The number of retries for the jobs. (Default: `3`)
- `forks`: How many jobs to run simultaneously? (Default: `1`)
- `scheduler`: The scheduler to use. If not specified, the scheduler specified in the pipeline level configurations will be used.
- `scheduler_opts`: The options for the scheduler. If not specified, the scheduler options specified in the pipeline level configurations will be used.
    - See also [`How to run the pipeline on a cluster?`](./faq.md#how-to-run-the-pipeline-on-a-cluster) for more details.

To know more about the configuration items for the pipeline, you can also read the [pipen docs](https://pwwang.github.io/pipen/configurations/).

## Enabling/disabling processes

By default, only the essential processes are enabled.

If scTCR-seq data is avaiable, these processes include:

- [`SampleInfo`](processes/SampleInfo.md)
- [`ImmunarchLoading`](processes/ImmunarchLoading.md)
- [`Immunarch`](processes/Immunarch.md)
- [`SeuratPreparing`](processes/SeuratPreparing.md)
- [`SeuratClustering`](processes/SeuratClustering.md)
- [`IntegratingTCR`](processes/IntegratingTCR.md)
- [`ClusterMarkers`](processes/ClusterMarkers.md)
- [`SeuratClusterStats`](processes/SeuratClusterStats.md)

If only scRNA-seq data is available, these processes include:

- [`SampleInfo`](processes/SampleInfo.md)
- [`SeuratPreparing`](processes/SeuratPreparing.md)
- [`SeuratClustering`](processes/SeuratClustering.md)
- [`ClusterMarkers`](processes/ClusterMarkers.md)
- [`SeuratClusterStats`](processes/SeuratClusterStats.md)

See also [`Routes of the pipeline`](./introduction.md#routes-of-the-pipeline) for more details.

To enable optional processes, you just need to add the corresponding sections for the processes in the configuration file. As long as the process name appears in the configuration file, the process will be enabled. For example, if you want to add module scores (e.g. cell activation score) to the `Seurat` object, you can add the following lines to the configuration file:

```toml
[ModuleScoreCalulator.envs.modules.TCell_Terminal_Differentiation]
features = ["TIGIT", "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "CD244", "CD160"]
```

If [`TCRClustering`](processes/TCRClustering.md) or [`TCRClusterStats`](processes/TCRClusterStats.md) is enabled, [`TCRClustering`](processes/TCRClustering.md), [`IntegratingTCRClusters`](processes/IntegratingTCRClusters.md), and [`TCRClusterStats`](processes/TCRClusterStats.md) will be enabled automatically.

/// Tip
You may find out that for some processes, the default configurations are good enough for you to run. For example, [`TCRClustering`](processes/TCRClustering.md) is not enabled by default. If you don't change any configurations (by not putting in the configuration file nor changing any items on the web interface of `pipen-board`) for the process, it will not be triggered. However, the default configurations are good enough for you to run the process. To enable it, you can either add this process manually in the configuration file:

```toml
# ... other configurations
[TCRClustering]
```

or if you are using `pipen-board`, you can change a configuration item that does not actually affect the process. For example, you can change the `forks` of the process to `2`, instead of the default `1`, since the process is a single-job process. Then the process will be put in the configuration file and will be enabled.

```toml
[TCRClustering]
forks = 2
```
///

## Minimal configurations

The minimal configurations are just the configurations with the input file:

```toml
[SampleInfo.in]
infile = [ "samples.txt" ]
```

The input file is the metadata file mentioned in [`Preparing the input`](./preparing-input.md#metadata).

With the minimal configurations, the pipeline will have [the essential processes](#enablingdisabling-processes) enabled, depending on whether scTCR-seq data is available or not.

You can also check the example report [here](http:/imp.pwwang.com/minimal/REPORTS/) to see what you will get with the minimal configurations, with scTCR-seq data available.

## Environment variable types

The types of environment variables are annotated in the brackets next the name of the environment variables. For example, the type of `envs.indicator_genes` of [`TCellSelection`](processes/TCellSelection.md) is `list`, and it's annotated as:

```
- indicator_genes (list): The genes to be used to select T cells.
```

By default, the type of environment variables is `string`. The annotated types are helpful for the environment variables to be passed from the command line. It defines the argument and helps parse the argument from the command line. It is also useful to define the input elements from the [`pipen-board`][1] web interface and parse the values passed from the web interface as desired types.

The following types are supported:

- `string`: The default type, the values will be used as strings.
- `int`: The values will be parsed as integers.
- `float`: The values will be parsed as floats.
- `flag`: The values will be parsed as boolean values.
- `list`/`array`: The values will be parsed as lists.
    - You can also see the `itype` of some environment variables, that specifies the type of the elements in the list. It must be atomatic types, such as `int`, `float`, `string`, and `flag`.
- `json`: The values will be reciived as JSON strings and parsed as dictionaries (in python).
- `choice`/`choices`: The value should be chosen from one of the choices listed as sub-items.
- `mchoice`/`mchoices`: The value should be chosen from one or more of the choices listed as sub-items.
- `ns`/`namespace`: There are sub-items for the value. The sub-items will be parsed as key-value pairs.

## Understanding the data

Understanding how the data is presented in the pipeline is helpful for the configuration, especially for the processes, such as [`RadarPlots`](processes/RadarPlots.md) and [`CellsDistribution`](processes/CellsDistribution.md). The configurations of this kind of processes are relying on the metadata.

You can refer to the individual process pages for more details. Here we just give an introduction of how it works to set the configurations.

### The assay of the `Seurat` object

The `Seurat` object is the main object used in the pipeline. You can have multiple assays in the `Seurat` object. While preparing the `Seurat` object at [`SeuratPreparing`](processes/SeuratPreparing.md) process, the default assay is determined. If `envs.use_sct` is true, meaning [`SCTransform`](https://satijalab.org/seurat/reference/sctransform) is used, the default assay will be `SCT`. If you are using `cca` or `rpca` integration, `integrated` will be used as the default assay. Otherwise, the default assay will be `RNA`. For downstream processes using the expression values, we provide an option to specify the assay to use. However, the default assay is used. Unless you know what you are doing, you are not recommended to change the default assay.

### Using existing columns in the metadata

In most cases, you can use the existing columns in the metadata to set the configurations. For example, if you want to plot the clone residency for each patient/subject, you need to specify the column name of the sample ID, as well as the column with the paired sample information (i.e. `tumor` vs `blood`).

Suppose the metadata (sitting in `immdata$meta` in `R` for example) is as follows:

|     Sample    | Source |
| ------------- | ------ |
| MM003-Eariler | BM     |
| MM003-Eariler | PB     |
| MM005-Eariler | BM     |
| MM005-Eariler | PB     |

Then you can set the configurations as follows:

```toml
[RadarPlots.envs]
subject = "Sample"
group = "Source"
```

And you will get the following plots:

![clone-residency](processes/images/CloneResidency.png)

### Mutating the metadata

Sometimes, you may want to mutate the metadata to get the desired information. Of course, you can have them prepared in the input file, as those extra columns with meta information will be attached to the object (either `immunarch$meta` or `srtobj@meta.data`) automatically. See [`Preparing the input`](./preparing-input.md) for more details. However, sometimes the metadata is specific to some processes, you may not want to have them prepared in the input file to get all processes contaminated. Moreover, those derived columns are usually based on the existing columns, so that is also helpful to create them on the fly to keep the input file clean.

In such case, for example, if you want to plot the clone residency for two groups (e.g. `BM-Pre` vs. `BM-Post`) of samples for the same group (e.g. `A`). However, the `Source` and `Timepoint` information are not in a single column of metadata. Here is when `mutaters` come in place.

Suppose the metadata (sitting in `immdata$meta` in `R` for example) is as follows:

| Sample | Group | Source | Timepoint |
| ------ | ----- | ------ | --------- |
| MM003  | A     | BM     | Pre       |
| MM003  | A     | BM     | Pre       |
| MM005  | A     | BM     | Post      |
| MM005  | A     | BM     | Post      |
| ...    | ...   | ...    | ...       |

Then you can set the configurations as follows:

```toml
[CloneResidency.envs.mutaters]
SampleGroup = "paste0(Sample, '-', Timepoint)"

[CloneResidency.envs]
subject = "Group"
group = "SampleGroup"
order = ["BM-Pre", "BM-Post"]
```

Then you will get a clone residency plot for group `A` with `BM-Pre` as x-axis and `BM-Post` as y-axis.

The key-value pairs of `mutaters` are passed to [`dplyr::mutate()`][11] function. The actual code to mutate the metadata is as follows:

```r
df %>% mutate(SampleGroup = paste0(Sample, '-', Timepoint))
```

So, for this kind of advanced configurations, you need to have some knowledge of `dplyr` in `R`.

You also need to pay attention to the keys of `mutaters`. Basically, the keys are the column names you want to create. So you need to make sure that the column names are not in the metadata already. Otherwise, the existing columns will be overwritten.

For scRNA-seq data, the existing column names of the metadata are:

- `orig.ident`
- `nCount_RNA`
- `nFeature_RNA`

and the meta columns in the input file. See also [`Preparing the input`](./preparing-input.md) for more details.

There could also be some other columns, depending on the previous processes. For example, if you have the cells clustered, there will be a column named `seurat_clusters` in the metadata.

For scTCR-seq data, `Sample` is the only existing column in the metadata after loaded. Then the meta columns from the input file will be attached to the metadata.

The best practice is to use a prefix for the column names you want to create. For example, if you want to create a column named `Sample`, you can use `my_Sample` instead. Then you can make sure that the column names are not in the metadata already.

The other thing you need to pay attention to is that you should try to avoid `.` or `-` in the column names. For example, if you want to create a column named `Sample-Source`, you can use `Sample_Source` instead. This is because that the column names will be used as the keys of the environment variables, and some processes will translate `-` into `.`. See also [`Namespace environment variables`](#namespace-environment-variables) for more details.

#### Mutater helpers

Other than the direct expressions that are parsed by [`rlang::parse_expr()`](https://rlang.r-lib.org/reference/parse_expr.html), for processes with `envs.mutaters`, we also provide some helper functions to make it easier to create the columns, especially for identifying the clones that are expanded, collapsed, vanished and emerged between two groups. The helper functions are:

- `expanded()`: Identify the expanded clones between two groups.
- `collapsed()`: Identify the collapsed clones between two groups.
- `vanished()`: Identify the vanished clones between two groups.
- `emerged()`: Identify the emerged clones between two groups.

The helper functions take the following arguments:

* `df`: The cell-level data. When used in `dplyr::mutate()`, you can use `.` to refer to the dataframe.
* `group.by`: The column name in metadata to group the cells.
* `idents`: The first group or both groups of cells to compare (value in `group-by` column). If only the first group is given, the rest of the cells (with non-NA in `group-by` column) will be used as the second group.
* `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).
* `each`: A column name (without quotes) in metadata to split the cells.
    Each comparison will be done for each value in this column.
* `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).
* `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.
    If numeric column is given, the values should be the same for all cells in the same group.
    This will not be checked (only the first value is used).
* `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.
* `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.
* `order`: The expression passed to `dplyr::arrange()` to order intermediate dataframe and get the ids in order accordingly.
    The intermediate dataframe includes the following columns:
    * `<id>`: The ids of clones (i.e. `CDR3.aa`).
    * `<each>`: The values in `each` column.
    * `ident_1`: The size of clones in the first group.
    * `ident_2`: The size of clones in the second group.
    * `.diff`: The difference between the sizes of clones in the first and second groups.
    * `.sum`: The sum of the sizes of clones in the first and second groups.
    * `.predicate`: Showing whether the clone is expanded/collapsed/emerged/vanished.
* `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.
* `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.

The returned values of the functions depend on the `debug` and `uniq` arguments. If `debug` is `TRUE`, the intermediate data frame will be returned. Otherwise, the ids will be returned.
If `uniq` is `TRUE`, the ids will be unique. Otherwise, the ids will be the same length as the number of rows in input data frame, which is useful for mutating the metadata.

Let's say we have the following data frame:

| CDR3.aa | Group |
| ------- | ----- |
| Clone1  | A     |
| Clone1  | A     |
| Clone2  | A     |
| Clone3  | A     |
| Clone3  | A     |
| Clone3  | A     |
| Clone4  | A     |
| Clone4  | A     |
| Clone2  | B     |
| Clone2  | B     |
| Clone3  | B     |
| Clone3  | B     |
| Clone4  | B     |
| Clone4  | B     |
| Clone4  | B     |

To identify the expanded clones between groups `A` and `B`, the intermediate data frame will be like:

```r
expanded(df, Group, "A", debug = TRUE)
# A tibble: 4 × 6
  CDR3.aa ident_1 ident_2 .predicate  .sum .diff
  <chr>     <int>   <int> <lgl>      <int> <int>
1 Clone3        3       2 TRUE           5     1
2 Clone4        2       3 FALSE          5    -1
3 Clone2        1       2 FALSE          3    -1
4 Clone1        2       0 FALSE          2     2
```

So the expanded clones are `Clone3`, and if you want to include the emerged clones, the result will be `Clone3` and `Clone1`.

```r
expanded(df, Group, "A")
[1] "Clone3"

expanded(df, Group, "A", include_emerged = TRUE)
[1] "Clone3" "Clone1"

# Change the order based on the difference
expanded(df, Group, "A", include_emerged = TRUE, order = desc(.diff))
[1] "Clone1" "Clone3"
```

If you want to add a column named `Expanded` to the metadata to save the clone ids:

```r
df %>% mutate(Expanded = expanded(df, Group, "A", uniq = FALSE))
# A tibble: 15 × 3
   CDR3.aa Group Expanded
   <chr>   <chr> <chr>
 1 Clone1  A     NA
 2 Clone1  A     NA
 3 Clone2  A     NA
 4 Clone3  A     Clone3
 5 Clone3  A     Clone3
 6 Clone3  A     Clone3
 7 Clone4  A     NA
 8 Clone4  A     NA
 9 Clone2  B     NA
10 Clone2  B     NA
11 Clone3  B     Clone3
12 Clone3  B     Clone3
13 Clone4  B     NA
14 Clone4  B     NA
15 Clone4  B     NA
```

### Filtering/Subsetting the data

In most processes where we need to filter the data, we don't provide an option for you to set the expression for [`dplyr::filter()`][14]. Instead, you can make use of the `mutaters` to create a column for filtering. For example, if you only want to plot clone residency for only one patient/subject (e.g. `MM003-Eariler`) in [`CloneResidency`](processes/CloneResidency.md), you can set the configurations as follows (suppose we have `Sample` and `Source` columns in the metadata):

```toml
[RadarPlots.envs.mutaters]
SingleSample = "if_else(Sample == 'MM003-Eariler', Sample, NA)"

[RadarPlots.envs]
subject = "SingleSample"
group = "Source"
```

Then you will get only one plot for `MM003-Eariler`, but not for `MM005-Eariler`. The `NA`s will be filtered out automatically.

## Namespace environment variables

There are some enviroment variables marked as `namespace`, which means that you can have sub-keys for them. For example, the `envs.SCTransform` of [`SeuratClusteringOfAllCells`](processes/SeuratClusteringOfAllCells.md) process is a `namespace` environment variable. It takes the arguments of [`Seurat::SCTransform()`][12] function. The names of arguments have dot (`.`) in them, such as `do.scale`, `do.center`, `seed.use`, etc. In the configuration file, we need to use dash (`-`) instead of dot (`.`) to set the values for these arguments. For example, if we want to set `do.scale` to `TRUE`, we need to set `do-scale` to `true` in the configuration file.

```toml
[SeuratClusteringOfAllCells.envs.SCTransform]
do-scale = true
```

This is because that we use [`pipen-args`][5] plugin backended by [`argx`][13] to parse the command line arguments, including the configuration file. If we use `.` directly in the configuration file:

```toml
[SeuratClusteringOfAllCells.envs.SCTransform]
do.scale = true
```

Then the `pipen-args` will parse it as `do` is the key and `scale` is the sub-key, and the above configuration will be parsed as:

```toml
[SeuratClusteringOfAllCells.envs.SCTransform]
do = { scale = true }
```

which is not what we want.

The reason why `.` is parsed as sub-key is that we want the argument to be able to be passed from command line. For example, if we want to set `do.scale` to `TRUE` from command line, we can do:

```bash
$ immunopipe --SeuratClusteringOfAllCells.envs.SCTransform.do-scale true
```

If we use `.` instead of `-`:

```bash
$ immunopipe --SeuratClusteringOfAllCells.envs.SCTransform.do.scale true
```

Then the `pipen-args` will parse it as

```toml
[SeuratClusteringOfAllCells.envs.SCTransform]
do = { scale = true }
```

again.

/// Tip
You don't need to worry about which environment variables are `namespace` ones. We will mention it in the individual process pages and the description of the environment variables in `pipen-board` configuration descriptions.
///

## Multi-case variable design

Some environment variables are designed to support multiple cases. However, in most cases, we only need to set the values for the default case. In such cases, the environment variable is usually a `namespace` environment variable with the sub-keys needed for the default case. In order to support multiple cases, a sub-key `cases` is added to the `namespace` environment variable. The `cases` is a dictionary (key-value pairs), where the keys are the names of the cases, and the values are the sub-keys for the corresponding cases. For example, the `envs.cluster_size` of [`TCRClusterStats`](processes/TCRClusterStats.md) process:

```toml
[TCRClusterStats.envs.cluster_size]
by = "Sample"
devpars = { width = 1000, height = 1000, res = 100 }
cases = {}
```

If `cases` is empty, the default case will be added automatically. The name of the default case is `DEFAULT`. So the above configuration is equivalent to:

```toml
[TCRClusterStats.envs.cluster_size]
by = "Sample"
devpars = { width = 1000, height = 1000, res = 100 }
cases = { DEFAULT = {} }
```

If you want to add more cases, you can add them to the `cases` dictionary. For example, if you want to add a case named `CASE1`, you can do:

```toml
[TCRClusterStats.envs.cluster_size]
by = "Sample"
devpars = { width = 1000, height = 1000, res = 100 }
cases = { DEFAULT = {}, CASE1 = {} }
```

Then you can set the values for the default case and `CASE1` case. For example, if you want to set the `by` column to `Sample` for the default case and `Sample1` for `CASE1`, you can do:

```toml
[TCRClusterStats.envs.cluster_size]
by = "Sample"
devpars = { width = 1000, height = 1000, res = 100 }
cases = { DEFAULT = { }, CASE1 = { by = "Sample1" } }
```

If a key in a case is not specified, the value in the default case will be used. In the above example,
`"Sample"` will be used for `by` of the `DEFAULT` case, and `{ width = 1000, height = 1000, res = 100 }` will be used for `devpars` of the `CASE1` case.

## Security alert

/// Danger
Note that some configuration items will be evaluated in the scripts directly. For example, the `mutaters` will be passed to `R` scripts, parsed and evaluated so that they can be used in `dplyr::mutate()`. Even though some were evaluated by [`rlang`][15], not all of them are safe. Some of them are evaluated directly. For example, one could inject malicious code in the expressions passed by `dplyr::filter()`. For example, in the script:

```r
df %>% filter({{ expression }})
```

The expected `expression` is something like `Sample == "Sample001"`. However, one could pass `Sample == "Sample001"); system("cat /etc/passwd")` to the `expression`, which will be evaluated as:

```r
df %>% filter(Sample == "Sample001"); system("cat /etc/passwd")
```

This will cause the pipeline to run the command `cat /etc/passwd` in the shell. This is just an example. One could do more harm by injecting malicious code.

When you give acess of composing the configuration file to others or the public (not recommended), either via the command line or the web interface by `pipen-board`, you need to be careful about the security issues.
///


[1]: https://github.com/pwwang/pipen-board
[2]: https://github.com/pwwang/pipen-verbose
[3]: https://github.com/pwwang/pipen-runinfo
[4]: https://github.com/pwwang/pipen-filters
[5]: https://github.com/pwwang/pipen-args
[6]: https://github.com/pwwang/pipen-annotate
[7]: https://github.com/pwwang/pipen-report
[8]: https://github.com/pwwang/pipen-log2file
[9]: https://github.com/pwwang/pipen-cli-run
[10]: https://github.com/pwwang/xqute
[11]: https://dplyr.tidyverse.org/reference/mutate.html
[12]: https://satijalab.org/seurat/reference/sctransform
[13]: https://github.com/pwwang/argx
[14]: https://dplyr.tidyverse.org/reference/filter.html
[15]: https://rlang.r-lib.org/
