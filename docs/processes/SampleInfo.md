# SampleInfo

List sample information and perform statistics



## Input

- `infile`:
    The input file to list sample information
    The input file should be a csv/tsv file with header

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
- `exclude_cols` *(`auto`)*: *Default: `TCRData,BCRData,RNAData`*. <br />
    The columns to exclude in the table in the report.<br />
    Could be a list or a string separated by comma.<br />
- `defaults` *(`ns`)*:
    The default parameters for `envs.stats`.<br />
    - `plot_type`: *Default: `bar`*. <br />
        The type of the plot.<br />
        See the supported plot types here:<br />
        <https://user.github.io/plotthis/reference/index.html>
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

