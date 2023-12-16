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
    - `distinct`:
- `stats` *(`type=json`)*: *Default: `{}`*. <br />
    The statistics to perform.<br />
    The keys are the case names and the values are the parameters
    inheirted from `envs.defaults`.<br />

## Examples

### Example data

| Subject | Sample | Source | Score |
| ------- | ------ | ------ | ----- |
| A       | A1     | Tumor  | 1     |
| A       | A2     | Numor  | 8     |
| A       | A3     | Tumor  |3      |
| A       | A4     | Normal |8      |
| B       | B1     | Tumor  |2      |
| B       | B2     | Normal |8      |
| B       | B3     | Tumor  |4      |
| B       | B4     | Normal |8      |
| C       | C1     | Tumor  |9      |
| C       | C2     | Normal |3      |
| C       | C3     | Tumor  |7      |
| C       | C4     | Normal |3      |
| D       | D1     | Tumor  |10     |
| D       | D2     | Normal |5      |
| D       | D3     | Tumor  |10     |
| D       | D4     | Normal |5      |
| E       | E1     | Tumor  |6      |
| E       | E2     | Normal |5      |
| E       | E3     | Tumor  |6      |
| E       | E4     | Normal |5      |
| F       | F1     | Tumor  |8      |
| F       | F2     | Normal |10     |
| F       | F3     | Tumor  |8      |
| F       | F4     | Normal |10     |

### Count the number of samples per Source

```toml
[SampleInfo.envs.stats]
Samples_Source = { "group": "Source" }
Samples_Source_each_Subject = { "group": "Source", "each": "Subject" }
```

![Samples_Source](../processes/images/SampleInfo_Samples_Source.png)

### Explore the distribution of the Score

```toml
[SampleInfo.envs.stats.Score_Source_vlnbox]
on = "Score"
group = "Source"
plot = "violin+box"

[SampleInfo.envs.stats.Score_Source_each_Subject_vlnbox]
on = "Score"
group = "Source"
plot = "violin+box"
each = "Subject"
```

![Score_Source](../processes/images/SampleInfo_Score_Source.png)

