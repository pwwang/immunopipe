# CellCellCommunicationPlots

Visualization for cell-cell communication inference.

R package [`CCPlotR`](https://github.com/Sarah145/CCPlotR) is used to visualize
the results.<br />

## Input

- `cccfile`:
    The output file from `CellCellCommunication`
    or a tab-separated file with the following columns: `source`, `target`,
    `ligand`, `receptor`, and `score`.<br />
    If so, `in.expfile` can be provided where `exp_df` is needed.<br />
- `expfile`:
    The expression file with the expression of ligands and receptors.<br />
    Columns include: `cell_type`, `gene` and `mean_exp`.<br />

## Output

- `outdir`: *Default: `{{in.cccfile | stem}}-ccc_plots`*. <br />
    The output directory for the plots.<br />

## Environment Variables

- `score_col`: *Default: `mag_score`*. <br />
    The column name in the input file that contains the score, if
    the input file is from `CellCellCommunication`.<br />
    Two alias columns are added in the result file of `CellCellCommunication`,
    `mag_score` and `spec_score`, which are the magnitude and specificity
    scores.<br />
- `subset`:
    An expression to pass to `dplyr::filter()` to subset the ccc data.<br />
- `cases` *(`type=json`)*: *Default: `{}`*. <br />
    The cases for the plots.<br />
    The keys are the names of the cases and the values are the arguments for
    the plots. The arguments include:<br />
    * kind: one of `arrow`, `circos`, `dotplot`, `heatmap`, `network`,
    and `sigmoid`.<br />
    * devpars: The parameters for `png()` for the plot, including `res`,
    `width`, and `height`.<br />
    * section: The section name for the report to group the plots.<br />
    * <other>: Other arguments for `cc_<kind>` function in `CCPlotR`.<br />
    See the documentation for more details.<br />
    Or you can use `?CCPlotR::cc_<kind>` in R.<br />

