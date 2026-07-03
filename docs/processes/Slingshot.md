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
- `dims` *(`type=auto`)*: *Default: `[1, 2]`*. <br />
    The dimensions to use for the analysis.<br />
    A list or a string with comma separated values.<br />
    Consecutive numbers can be specified with a colon (`:`) or a dash (`-`).<br />
- `start`:
    The starting group for the Slingshot analysis.<br />
- `end`:
    The ending group for the Slingshot analysis.<br />
- `prefix`:
    The prefix to add to the column names of the resulting pseudotime variable.<br />
- `reverse` *(`flag`)*: *Default: `False`*. <br />
    Logical value indicating whether to reverse the pseudotime variable.<br />
- `align_start` *(`flag`)*: *Default: `False`*. <br />
    Whether to align the starting pseudotime values at the maximum pseudotime.<br />
- `seed` *(`type=int`)*: *Default: `8525`*. <br />
    The seed for the random number generator.<br />

