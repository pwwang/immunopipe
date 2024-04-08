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

![infile](../processes/images/SampleInfo-infile.png)

Theroetically, we can have multiple input files. However, it is not tested yet.<br />
If you have multiple input files to run, please run it with a different pipeline
instance (configuration file).<br />

For the content of the input file, please see details
[here](../preparing-input.md#metadata).<br />

Once the pipeline is finished, you can see the sample information in the report

![report](../processes/images/SampleInfo-report.png)

Note that the required `RNAData` and `TCRData` columns are not shown in the report.<br />
They are used to specify the paths of the `scRNA-seq` and `scTCR-seq` data, respectively.<br />

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
- `exclude_cols`: *Default: `TCRData,RNAData`*. <br />
    The columns to exclude in the table in the report.<br />
    Could be a list or a string separated by comma.<br />
- `defaults` *(`ns`)*:
    The default parameters for `envs.stats`.<br />
    - `on`: *Default: `Sample`*. <br />
        The column name in the data for the stats.<br />
        Default is `Sample`. The column could be either continuous or not.<br />
    - `subset`:
        An R expression to subset the data.<br />
        If you want to keep the distinct records, you can use
        `!duplicated(<col>)`.<br />
    - `group`:
        The column name in the data for the group ids.<br />
        If not provided, all records will be regarded as one group.<br />
    - `na_group` *(`flag`)*: *Default: `False`*. <br />
        Whether to include `NA`s in the group.<br />
    - `each`:
        The column in the data to split the analysis in different
        plots.<br />
    - `ncol` *(`type=int`)*: *Default: `2`*. <br />
        The number of columns in the plot when `each`
        is not `NULL`. Default is 2.<br />
    - `na_each` *(`flag`)*: *Default: `False`*. <br />
        Whether to include `NA`s in the `each` column.<br />
    - `plot`:
        Type of plot. If `on` is continuous, it could be
        `boxplot` (default), `violin`, `violin+boxplot` or `histogram`.<br />
        If `on` is not continuous, it could be `barplot` or
        `pie` (default).<br />
    - `devpars` *(`ns`)*:
        The device parameters for the plot.<br />
        - `width` *(`type=int`)*: *Default: `800`*. <br />
            The width of the plot.<br />
        - `height` *(`type=int`)*: *Default: `600`*. <br />
            The height of the plot.<br />
        - `res` *(`type=int`)*: *Default: `100`*. <br />
            The resolution of the plot.<br />
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
on = "sample"
group = "Diagnosis"
```

![Samples_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/dev/tests/output/sampleinfo/SampleInfo/N_Samples_per_Diagnosis (pie).png)

What if we want a bar plot instead of a pie chart?<br />

```toml
[SampleInfo.envs.stats."N_Samples_per_Diagnosis (bar)"]
on = "sample"
group = "Diagnosis"
plot = "barplot"
```

![Samples_Diagnosis_bar](https://raw.githubusercontent.com/pwwang/immunopipe/dev/tests/output/sampleinfo/SampleInfo/N_Samples_per_Diagnosis (bar).png)

### Explore Age distribution

The distribution of Age of all samples

```toml
[SampleInfo.envs.stats."Age_distribution (boxplot)"]
on = "Age"
```

![Age_distribution](https://raw.githubusercontent.com/pwwang/immunopipe/dev/tests/output/sampleinfo/SampleInfo/Age_distribution (boxplot).png)

How about the distribution of Age in each Diagnosis, and make it violin + boxplot?<br />

```toml
[SampleInfo.envs.stats."Age_distribution_per_Diagnosis (violin + boxplot)"]
on = "Age"
group = "Diagnosis"
plot = "violin+boxplot"
```

![Age_distribution_per_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/dev/tests/output/sampleinfo/SampleInfo/Age_distribution_per_Diagnosis (violin + boxplot).png)

How about Age distribution per Sex in each Diagnosis?<br />

```toml
[SampleInfo.envs.stats."Age_distribution_per_Sex_in_each_Diagnosis (boxplot)"]
on = "Age"
group = "Sex"
each = "Diagnosis"
plot = "boxplot"
ncol = 3
devpars = {height = 450}
```

![Age_distribution_per_Sex_in_each_Diagnosis](https://raw.githubusercontent.com/pwwang/immunopipe/dev/tests/output/sampleinfo/SampleInfo/Age_distribution_per_Sex_in_each_Diagnosis (boxplot).png)

