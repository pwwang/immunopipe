# CellCellCommunicationPlots

Visualization for cell-cell communication inference.



## Input

- `cccfile`:
    The output file from `CellCellCommunication`

## Output

- `outdir`: *Default: `{{in.cccfile | stem}}_plots`*. <br />
    The output directory for the plots.<br />

## Environment Variables

- `subset`:
    An expression to pass to `dplyr::filter()` to subset the ccc data.<br />
- `magnitude`:
    The column name in the data to use as the magnitude of the
    communication. By default, the second last column will be used.<br />
    See `li.mt.show_methods()` for the available methods in LIANA. or
    <https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html#Tileplot>
- `specificity`:
    The column name in the data to use as the specificity of the communication.<br />
    By default, the last column will be used. If the method doesn't have a specificity, set it to None.<br />
- `devpars` *(`ns`)*:
    The parameters for the plot.<br />
    - `res` *(`type=int`)*: *Default: `100`*. <br />
        The resolution of the plot
    - `height` *(`type=int`)*:
        The height of the plot
    - `width` *(`type=int`)*:
        The width of the plot
- `more_formats` *(`type=list`)*: *Default: `[]`*. <br />
    The additional formats to save the plots.<br />
- `descr`: *Default: `Cell-cell communication plot`*. <br />
    The description of the plot.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    The cases for the plots.<br />
    The keys are the names of the cases and the values are the arguments for
    the plots. The arguments include the ones inherited from `envs`.<br />
- `<more>`:
    Other arguments passed to
    [scplotter::CCCPlot](https://pwwang.github.io/scplotter/reference/CCCPlot.html)

