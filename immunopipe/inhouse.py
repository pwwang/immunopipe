"""In house processes"""
from biopipen.core.proc import Proc
from biopipen.core.config import config
from pipen_annotate import annotate

TEST_OUTPUT_BASEURL = "https://raw.githubusercontent.com/pwwang/immunopipe/tests-output"


@annotate.format_doc(vars={"output_baseurl": TEST_OUTPUT_BASEURL})
class TOrBCellSelection(Proc):
    """Separate T and non-T cells and select T cells; or
    separate B and non-B cells and select B cells.

    If all of your cells are T/B cells, do not set any configurations for this process.

    In such a case, [`SeuratClusteringOfAllCells`](SeuratClusteringOfAllCells.md) should
    not be used, and [`SeuratClustering`](SeuratClustering.md) will be clustering all
    of the cells, which are all T/B cells.

    There are two ways to separate T and non-T cells; or B and non-B cells:

    1. Use the an expression indicator directly from the metadata.
    2. Use the expression values of indicator genes, and the clonotype percentage
    of the clusters.

    You can also use indicator gene expression values only to select T/B cells by
    setting `envs.ignore_vdj` to true.

    Examples:

        ### Use T cell indicator directly

        If you have a metadata like this:

        | id | Clonotype_Pct | seurat_clusters |
        |----|---------------|-----------------|
        | 1  | 0.1           | 1               |
        | 2  | 0.3           | 2               |
        | 3  | 0.5           | 3               |

        With the configuration below:

        ```toml
        [TOrBCellSelection.envs]
        selector = "Clonotype_Pct > 0.25"
        ```

        The T cells will be selected as:

        | id | Clonotype_Pct | seurat_clusters | is_TCell |
        |----|---------------|-----------------|----------|
        | 1  | 0.1           | 1               | FALSE    |
        | 2  | 0.3           | 2               | TRUE     |
        | 3  | 0.5           | 3               | TRUE     |

        ### Use indicator genes

        Let's say we set the indicator genes to `["CD3D", "CD3E", "CD3G"]`.

        The mean expression values will be calculated for each cluster:

        | id | Clonotype_Pct | seurat_clusters | CD3D | CD3E | CD3G |
        |----|---------------|-----------------|------|------|------|
        | 1  | 0.1           | 1               | 0.1  | 0.0  | 0.1  |
        | 2  | 0.3           | 2               | 1.2  | 1.3  | 0.6  |
        | 3  | 0.5           | 3               | 1.5  | 0.8  | 0.9  |

        Then a kmeans clustering will be performed on the mean expression values of
        the indicator genes, together with `Clonotype_Pct`, with K=2.

        | id | Clonotype_Pct | seurat_clusters | CD3D | CD3E | CD3G | is_TCell |
        |----|---------------|-----------------|------|------|------|----------|
        | 1  | 0.1           | 1               | 0.1  | 0.0  | 0.1  | FALSE    |
        | 2  | 0.3           | 2               | 1.2  | 1.3  | 0.6  | TRUE     |
        | 3  | 0.5           | 3               | 1.5  | 0.8  | 0.9  | TRUE     |

        ![kmeans](images/TOrBCellSelection-kmeans.png)

        The cluster with higher clonoype percentage will be selected as T/B cells
        (`is_selected = TRUE`), and sent to
        [`SeuratClustering`](SeuratClustering.md) for
        further clustering and downstream analysis.

        ### Details of the selection

        More details are saved in the `details` output directory, for example
        the indicator gene expression plots and the selected cells:

        ![feature plots]({{output_baseurl}}/torbcellselection/TOrBCellSelection/details/feature-plots.png){: width="80%" }
        ![selected cells]({{output_baseurl}}/torbcellselection/TOrBCellSelection/details/selected-cells.png){: width="80%" }
        ![selected cell fractions]({{output_baseurl}}/torbcellselection/TOrBCellSelection/details/selected_cells_pie.png){: width="80%" }
        ![CD3D vs CD3E]({{output_baseurl}}/torbcellselection/TOrBCellSelection/details/CD3D-vs-CD3E.png){: width="80%" }
        ![CD3G vs CD3E]({{output_baseurl}}/torbcellselection/TOrBCellSelection/details/CD3G-vs-CD3E.png){: width="80%" }

        The same plots are available when the cells are selected from a file of
        cell barcodes:

        ![feature plots (cell barcodes)]({{output_baseurl}}/torbcellselection_with_cellbarcodes/TOrBCellSelection/details/feature-plots.png){: width="80%" }
        ![selected cells (cell barcodes)]({{output_baseurl}}/torbcellselection_with_cellbarcodes/TOrBCellSelection/details/selected_cells_per_sample.png){: width="80%" }
        ![selected cell fractions (cell barcodes)]({{output_baseurl}}/torbcellselection_with_cellbarcodes/TOrBCellSelection/details/selected_cells_pie.png){: width="80%" }

    Input:
        srtobj: Seurat object file in RDS/qs2
        immdata: Immune repertoire data file in RDS/qs2

    Output:
        outfile: Seurat object file in qs2 format
        outdir: Output directory with details

    Envs:
        ignore_vdj (flag): Ignore VDJ information for T/B cell selection.
            Use only the expression values of indicator genes if True.
            In this case, the `Clonotype_Pct` column does not exist in the metadata.
            If you want to use `k-means` to select T/B cells, you must have more than
            1 indicator gene, and the first indicator gene in `envs.indicator_genes`
            must be a positive marker, which will be used to select the cluster with
            higher expression values as T/B cells.
        selector: The expression passed to `tidyseurat::mutate(is_TCell = ...)`
            to indicate whether a cell is a T cell, or
            a file containing the cell barcodes regarding as T cells
            (beginning with `file://`).
            For example, `Clonotype_Pct > 0.25`
            to indicate cells with clonotype percentage > 25% are T cells.
            If `indicator_genes` is provided, the expression values can also be used
            in the expression. For example, `Clonotype_Pct > 0.25 & CD3E > 0`.
            If `selector` is not provided, a kmeans clustering will be performed
            on the expression values of `indicator_genes` and `Clonotype_Pct`,
            with K=2, and the cluster with higher clonotype percentage will be selected
            as T/B cells.
            With a file is provided, the column name of the cell barcodes can be
            provided after the path to the file, separated by `#`. For example,
            `file://path/to/cell_barcodes.txt#cell_barcode` to indicate
            the cell barcodes are in the column `cell_barcode` of the file.
            If the column name is not provided, the file is assumed to have only
            one column without a header, and the cell barcodes will be read from
            the first column.
        indicator_genes (list): A list of indicator genes whose expression values and
            clonotype percentage will be used to determine T/B cells.
            The markers could be either positive, such as `CD3E`, `CD3D`, `CD3G`, or
            negative, such as `CD19`, `CD14`, `CD68`, for T cells. For B cells,
            markers such as `CD19`, `MS4A1` (CD20), `CD79A`, `CD79B` could be used.

        kmeans (type=json): The parameters for `kmeans` clustering.
            Other arguments for [`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html)
            can be provided here. If there are dots in the argument names, replace them
            with `-`.
    """  # noqa: E501
    input = "srtobj:file, immdata:file"
    output = "outfile:file:{{in.srtobj | stem}}.qs, outdir:dir:details"
    envs = {
        "ignore_vdj": False,
        "selector": None,
        "indicator_genes": ["CD3E"],
        "kmeans": {"nstart": 25},
    }
    lang = config.lang.rscript
    script = "file://scripts/TOrBCellSelection.R"
    plugin_opts = {"report": "file://reports/TOrBCellSelection.svelte"}
    order = 5
