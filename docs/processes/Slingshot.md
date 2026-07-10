# Slingshot

Trajectory inference using Slingshot

This process is implemented based on the R package `slingshot`.<br />

## Input

- `sobjfile`:
    The seurat object file in RDS or qs format.<br />

## Output

- `outfile`: *Default: `{{in.sobjfile | stem}}.qs`*. <br />
    The output object with the trajectory information.<br />
    The lineages are stored in the metadata of the seurat object at
    columns `LineageX`, where X is the lineage number. The `BranchID`
    column contains the branch id for each cell.<br />
    One can use
    `scplotter::CellDimPlot(object, lineages = c("Lineage1", "Lineage2", ...))`
    to visualize the trajectories.<br />

## Environment Variables

- `group_by`:
    The column name in metadata to group the cells.<br />
    Typically, this column should be the cluster id.<br />
    Default is the default identity of the seurat object.<br />
- `reduction`:
    The nonlinear reduction to use for the trajectory analysis.<br />
- `dims` *(`type=auto`)*:
    The dimensions to use for the analysis.<br />
    A list or a string with comma separated values.<br />
    Consecutive numbers can be specified with a colon (`:`) or a dash (`-`).<br />
    Or a single number greater than 1, which will be expanded to `1:number`.<br />
    If `None`, all dimensions will be used.<br />
- `start`:
    The starting group for the Slingshot analysis.<br />
- `end`:
    The ending group for the Slingshot analysis.<br />
- `reverse` *(`flag`)*: *Default: `False`*. <br />
    Logical value indicating whether to reverse the pseudotime variable.<br />
- `align_start` *(`flag`)*: *Default: `False`*. <br />
    Whether to align the starting pseudotime values at the maximum pseudotime.<br />
- `seed` *(`type=int`)*: *Default: `8525`*. <br />
    The seed for the random number generator.<br />
- `cases`: *Default: `{}`*. <br />
    A dictionary of cases to run the analysis.<br />
    The keys are the names of the cases, which will be served as the
    prefix to add to the column names of the resulting pseudotime variable.<br />
    For example, if the case name is `case1`, the resulting pseudotime variable
    will be stored in the column `case1_LineageX` and `case1_BranchID`.<br />
    The values are the arguments and will be inherited from the `envs` above, except for `cases`.<br />
    The default case will be added with the default values under `envs` with an empty prefix.<br />
- `outtype`: *Default: `qs2`*. <br />

