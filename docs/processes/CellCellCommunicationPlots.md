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
    You can have a special `plot_type` `"table"` to generate a table for the
    ccc data to save as a text file and show in the report.<br />
    If no cases are given, a default case will be used, with the
    key `Cell-Cell Communication`.<br />
- `<more>`:
    Other arguments passed to
    [scplotter::CCCPlot](https://pwwang.github.io/scplotter/reference/CCCPlot.html)

## Examples

### Network Plot

```toml
[CellCellCommunicationPlots.envs.cases."Cell-Cell Communication Network"]
plot_type = "network"
legend-position = "none"
theme = "theme_blank"
theme_args = {add_coord = false}
```

![Network Plot](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/cccplots/CellCellCommunicationPlots/sampleinfo.scRep-ccc_plots/Cell-Cell-Communication-Network.png){: width="80%"}

### Circos Plot

![Circos Plot](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/cccplots/CellCellCommunicationPlots/sampleinfo.scRep-ccc_plots/Cell-Cell-Communication-Circos-Plot.png){: width="80%"}

### Heatmap Plot

![Heatmap Plot](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/cccplots/CellCellCommunicationPlots/sampleinfo.scRep-ccc_plots/Cell-Cell-Communication-Heatmap.png){: width="80%"}

### Cell-Cell Communication Interaction (Box Plot)

```toml
[CellCellCommunicationPlots.envs.cases."Cell-Cell Communication Interaction (Box Plot)"]
plot_type = "box"
x_text_angle = 90
method = "interaction"
```

![Box Plot](https://raw.githubusercontent.com/pwwang/immunopipe/tests-output/cccplots/CellCellCommunicationPlots/sampleinfo.scRep-ccc_plots/Cell-Cell-Communication-Interaction-Box-Plot-.png){: width="80%"}

