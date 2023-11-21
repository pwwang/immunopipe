# SeuratMetadataMutater

Attach TCR clone information as meta columns to Seurat object

This process is used to integrate scTCR-seq data into the `Seurat` object.<br />
The scTCR-seq data is loaded by [ImmunarchLoading](./ImmunarchLoading.md) process.<br />
The integration is done by matching the barcodes from the `Seurat` object and
the scTCR-seq data.<br />
The barcodes from the scTCR-seq data are prefixed with the sample name,
for example, `Sample1_AAACCTGAGAAGGCTA-1`. The prefix is specified by the `prefix`
environment variable in the [ImmunarchLoading](./ImmunarchLoading.md) process.<br />

[ImmunarchLoading](./ImmunarchLoading.md) process will generate a text file with
the information for each cell.<br />
`ImmunarchLoading.envs.metacols` can be used to specify the columns to be exported
to the text file, which will then be integrated into the `Seurat` object
by this process.<br />

You may also use `envs.mutaters` to add new columns to the metadata.<br />
These columns can be used for downstream analysis.<br />
An additional column `TCR_Presence` is added so later on we can overlay the
TCR presence on the UMAP plot in
[`SeuratClusteringOfTCells`](./SeuratClusteringOfTCells.md) process.<br />

/// Warning
If you are modifying `envs.mutaters`, make sure you keep the `TCR_Presence` column.<br />
Because by default, `SeuratClusteringOfTCells` process will use this column to
overlay the TCR presence on the UMAP plot.<br />
///


## Environment Variables

- `mutaters` *(`type=json`)*: *Default: `{'TCR_Presence': 'if_else(is.na(CDR3.aa), "TCR_absent", "TCR_present")'}`*. <br />
    The mutaters to mutate the metadata.<br />
    The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.<br />
    There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`,
    which can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).<br />
    For example, you can use
    `{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq = FALSE)"}`
    to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
    with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1.<br />
    The values in this columns for other clones will be `NA`.<br />
    Those functions take following arguments:<br />
    * `df`: The metadata data frame. You can use the `.` to refer to it.<br />
    * `group-by`: The column name in metadata to group the cells.<br />
    * `idents`: The first group or both groups of cells to compare (value in `group-by` column). If only the first group is given, the rest of the cells (with non-NA in `group-by` column) will be used as the second group.<br />
    * `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).<br />
    * `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).<br />
    * `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.<br />
    If numeric column is given, the values should be the same for all cells in the same group.<br />
    This will not be checked (only the first value is used).<br />
    * `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.<br />
    * `order`: The order of the returned ids. It could be `sum` or `diff`, which is the sum or diff of the `compare` between idents.<br />
    Two kinds of modifiers can be added, including `desc` and `abs`.<br />
    For example, `sum,desc` means the sum of `compare` between idents in descending order.<br />
    Default is `diff,abs,desc`. It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
    ids will be in the same order as in `df`.<br />
    * `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.<br />
    * `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.<br />

