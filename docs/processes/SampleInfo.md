# SampleInfo

List sample information and perform statistics



This process is the entrance of the pipeline. It just pass by input file and list
the sample information in the report.<br />

To specify the input file in the configuration file, use the following

```toml
[SampleInfo.in]
infile = [ "path/to/sample_info.txt" ]
```

Or with `pipen-board`, find the `SampleInfo` process and click the `Edit` button.<br />
Then you can specify the input file here

![infile](../..//processes/images/SampleInfo-infile.png)

Multiple input files are supported by the underlying pipeline framework. However,
we recommend to run it with a different pipeline instance with configuration files.<br />

For the content of the input file, please see details
[here](../..//preparing-input.md#metadata).<br />

You can add some columns to the input file while doing the statistics or you can even
pass them on to the next processes. See `envs.mutaters` and `envs.save_mutated`.<br />
But if you are adding a factor (categorical) column with desired levels, the order
can't be guaranteed, because we are saving them to a text file, where we can't guarantee
the order of the levels. If you want to add a factor column with desired levels, you can
set `envs.mutaters` of the `SeuratPreparing` process to mutate the column.<br />

Once the pipeline is finished, you can see the sample information in the report

![report](../..//processes/images/SampleInfo-report.png)

Note that the required `RNAData` (if not loaded from a Seurat object) and `TCRData`
columns are not shown in the report.<br />
They are used to specify the paths of the `scRNA-seq` and `scTCR-seq` data, respectively.<br />
Also note that when `RNAData` is loaded from a Seurat object (specified in the
`LoadingRNAFromSeurat` process), the metadata provided in this process will not be
integrated into the Seurat object in the downstream processes. To incoporate
these meta information into the Seurat object, please provide them in the
Seurat object itself or use the `envs.mutaters` of the `SeuratPreparing` process
to mutate the metadata of the Seurat object. But the meta information provided in this
process can still be used in the statistics and plots in the report.<br />

You may also perform some statistics on the sample information, for example,
number of samples per group. See next section for details.<br />

/// Tip
This is the start process of the pipeline. Once you change the parameters for
this process, the whole pipeline will be re-run.<br />

If you just want to change the parameters for the statistics, and use the
cached (previous) results for other processes, you can set `cache` at
pipeline level to `"force"` to force the pipeline to use the cached results
and `cache` of `SampleInfo` to `false` to force the pipeline to re-run the
`SampleInfo` process only.<br />

```toml
cache = "force"

[SampleInfo]
cache = false
```
///

## Input

- `infile` *(`required`)*:
    The input file to list sample information The input file should be a csv/tsv file with header.<br />
    The input file should have the following columns.<br />
    * Sample: A unique id for each sample.<br />
    * TCRData: The directory for single-cell TCR data for this sample.<br />
    Specifically, it should contain filtered_contig_annotations.csv
    or all_contig_annotations.csv from cellranger.<br />
    * RNAData: The directory for single-cell RNA data for this sample.<br />
    Specifically, it should be able to be read by
    `Seurat::Read10X()` or `Seurat::Read10X_h5()` or `SeuratDisk::LoadLoom()`.<br />
    See also https://satijalab.org/seurat/reference/read10x.<br />
    * Other columns are optional and will be treated as metadata for
    each sample.<br />

## Output

- `outfile`: *Default: `{{in.infile | basename}}`*. <br />
    The output file with sample information, with mutated columns
    if `envs.save_mutated` is True.<br />
    The basename of the output file will be the same as the input file.<br />
    The file name of each plot will be slugified from the case name.<br />
    Each plot has 3 formats: pdf, png and code.zip, which contains the
    data and R code to reproduce the plot.<br />

## Environment Variables

- `sep`: *Default: `	`*. <br />
    The separator of the input file.<br />
- `mutaters` *(`type=json`)*: *Default: `{}`*. <br />
    A dict of mutaters to mutate the data frame.<br />
    The key is the column name and the value is the R expression
    to mutate the column. The dict will be transformed to a list in R
    and passed to `dplyr::mutate`.<br />
    You may also use `paired()` to identify paired samples. The function
    takes following arguments:<br />
    * `df`: The data frame. Use `.` if the function is called in
    a dplyr pipe.<br />
    * `id_col`: The column name in `df` for the ids to be returned in
    the final output.<br />
    * `compare_col`: The column name in `df` to compare the values for
    each id in `id_col`.<br />
    * `idents`: The values in `compare_col` to compare. It could be
    either an an integer or a vector. If it is an integer,
    the number of values in `compare_col` must be the same as
    the integer for the `id` to be regarded as paired. If it is
    a vector, the values in `compare_col` must be the same
    as the values in `idents` for the `id` to be regarded as paired.<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`.<br />
    If `FALSE`, you can mutate the meta data frame with the
    returned ids. Non-paired ids will be `NA`.<br />
- `save_mutated` *(`flag`)*: *Default: `False`*. <br />
    Whether to save the mutated columns.<br />
- `exclude_cols` *(`auto`)*: *Default: `TCRData,RNAData`*. <br />
    The columns to exclude in the table in the report.<br />
    Could be a list or a string separated by comma.<br />
- `defaults` *(`ns`)*:
    The default parameters for `envs.stats`.<br />
    - `plot_type`: *Default: `bar`*. <br />
        The type of the plot.<br />
        See the supported plot types here:<br />
        <https://pwwang.github.io/plotthis/reference/index.html>
        The plot_type should be lower case and the plot function used in
        `plotthis` should be used. The mapping from plot_type to the
        plot function is like `bar -> BarPlot`, `box -> BoxPlot`, etc.<br />
    - `more_formats` *(`list`)*: *Default: `[]`*. <br />
        The additional formats to save the plot.<br />
        By default, the plot will be saved in png, which is also used to
        display in the report. You can add more formats to save the plot.<br />
        For example, `more_formats = ["pdf", "svg"]`.<br />
    - `save_code` *(`flag`)*: *Default: `False`*. <br />
        Whether to save the R code to reproduce the plot.<br />
        The data used to plot will also be saved.<br />
    - `subset`:
        An expression to subset the data frame before plotting.<br />
        The expression should be a string of R expression that will be passed
        to `dplyr::filter`. For example, `subset = "Sample == 'A'"`.<br />
    - `section`:
        The section name in the report.<br />
        In case you want to group the plots in the report.<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plot.<br />
        - `width` *(`type=int`)*:
            The width of the plot.<br />
        - `height` *(`type=int`)*:
            The height of the plot.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plot.<br />
    - `descr`:
        The description of the plot, shown in the report.<br />
    - `<more>`:
        You can add more parameters to the defaults.<br />
        These parameters will be expanded to the `envs.stats` for each case,
        and passed to individual plot functions.<br />
- `stats` *(`type=json`)*: *Default: `{}`*. <br />
    The statistics to perform.<br />
    The keys are the case names and the values are the parameters
    inheirted from `envs.defaults`.<br />

## Examples

### Example data

| Sample | Age | Sex | Diagnosis |
|--------|-----|-----|-----------|
| C1     | 62  | F   | Colitis   |
| C2     | 71.2| F   | Colitis   |
| C3     | 56.2| M   | Colitis   |
| C4     | 61.5| M   | Colitis   |
| C5     | 72.8| M   | Colitis   |
| C6     | 78.4| M   | Colitis   |
| C7     | 61.6| F   | Colitis   |
| C8     | 49.5| F   | Colitis   |
| NC1    | 43.6| M   | NoColitis |
| NC2    | 68.1| M   | NoColitis |
| NC3    | 70.5| F   | NoColitis |
| NC4    | 63.7| M   | NoColitis |
| NC5    | 58.5| M   | NoColitis |
| NC6    | 49.3| F   | NoColitis |
| CT1    | 21.4| F   | Control   |
| CT2    | 61.7| M   | Control   |
| CT3    | 50.5| M   | Control   |
| CT4    | 43.4| M   | Control   |
| CT5    | 70.6| F   | Control   |
| CT6    | 44.3| M   | Control   |
| CT7    | 50.2| M   | Control   |
| CT8    | 61.5| F   | Control   |

### Count the number of samples per Diagnosis

```toml
[SampleInfo.envs.stats."N_Samples_per_Diagnosis (pie)"]
plot_type = "pie"
x = "sample"
split_by = "Diagnosis"
```

![Samples_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/sampleinfo/SampleInfo/N_Samples_per_Diagnosis-pie-.png)

What if we want a bar plot instead of a pie chart?<br />

```toml
[SampleInfo.envs.stats."N_Samples_per_Diagnosis (bar)"]
plot_type = "bar"
x = "Sample"
split_by = "Diagnosis"
```

![Samples_Diagnosis_bar](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/sampleinfo/SampleInfo/N_Samples_per_Diagnosis-bar-.png)

### Explore Age distribution

The distribution of Age of all samples

```toml
[SampleInfo.envs.stats."Age_distribution (histogram)"]
plot_type = "histogram"
x = "Age"
```

![Age_distribution](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/sampleinfo/SampleInfo/Age_distribution-Histogram-.png)

How about the distribution of Age in each Diagnosis, and make it violin + boxplot?<br />

```toml
[SampleInfo.envs.stats."Age_distribution_per_Diagnosis (violin + boxplot)"]
y = "Age"
x = "Diagnosis"
plot_type = "violin"
add_box = true
```

![Age_distribution_per_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/sampleinfo/SampleInfo/Age_distribution_per_Diagnosis-violin-boxplot-.png)

How about Age distribution per Sex in each Diagnosis?<br />

```toml
[SampleInfo.envs.stats."Age_distribution_per_Sex_in_each_Diagnosis (boxplot)"]
y = "Age"
x = "Sex"
split_by = "Diagnosis"
plot_type = "box"
ncol = 3
devpars = {height = 450}
```

![Age_distribution_per_Sex_in_each_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/sampleinfo/SampleInfo/Age_distribution_per_Sex_in_each_Diagnosis-boxplot-.png)

