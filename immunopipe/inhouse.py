"""In house processes"""
from biopipen.core.proc import Proc
from biopipen.core.config import config


class TCellSelection(Proc):
    """Separate T and non-T cells and select T cells.

    If all of your cells are T cells, do not set any configurations for this process.

    In such a case, [`SeuratClusteringOfAllCells`](SeuratClusteringOfAllCells.md) should
    not be used, and [`SeuratClustering`](SeuratClustering.md) will be clustering all
    of the cells, which are all T cells.

    There are two ways to separate T and non-T cells:

    1. Use the an expression indicator directly from the metadata.
    2. Use the expression values of indicator genes, and the clonotype percentage
    of the clusters.

    You can also use indicator gene expression values only to select T cells by setting
    `envs.ignore_tcr` to true.

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
        [TCellSelection.envs]
        tcell_selector = "Clonotype_Pct > 0.25"
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

        ![kmeans](images/TCellSelection-kmeans.png)

        The cluster with higher clonoype percentage will be selected as T cells
        (`is_TCell = TRUE`), and sent to
        [`SeuratClustering`](SeuratClustering.md) for
        further clustering and downstream analysis.

    Input:
        srtobj: Seurat object file in RDS
        immdata: Immunarch data file in RDS

    Output:
        rdsfile: Seurat object file in RDS
        outdir: Output directory with details

    Envs:
        ignore_tcr (flag): Ignore TCR information for T cell selection.
            Use only the expression values of indicator genes.
            In this case, the `Clonotype_Pct` column does not exist in the metadata.
            If you want to use `k-means` to select T cells, you must have more than
            1 indicator gene, and the first indicator gene in `envs.indicator_genes`
            must be a positive marker, which will be used to select the cluster with
            higher expression values as T cells.
        tcell_selector: The expression passed to `tidyseurat::mutate(is_TCell = ...)`
            to indicate whether a cell is a T cell. For example, `Clonotype_Pct > 0.25`
            to indicate cells with clonotype percentage > 25% are T cells.
            If `indicator_genes` is provided, the expression values can also be used
            in the expression. For example, `Clonotype_Pct > 0.25 & CD3E > 0`.
            If not provided, a kmeans clustering will be performed on the expression
            values of `indicator_genes` and `Clonotype_Pct`, with K=2, and the cluster
            with higher clonotype percentage will be selected as T cells.

            /// Tip | Changed in `0.11.0`
            `envs.tcell_indicator` is renamed to `envs.tcell_selector` in `0.11.0`.
            ///

        indicator_genes (list): A list of indicator genes whose expression values and
            clonotype percentage will be used to determine T cells.
            The markers could be either positive, such as `CD3E`, `CD3D`, `CD3G`, or
            negative, such as `CD19`, `CD14`, `CD68`.

        kmeans (type=json): The parameters for `kmeans` clustering.
            Other arguments for [`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html)
            can be provided here. If there are dots in the argument names, replace them
            with `-`.
    """
    input = "srtobj:file, immdata:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS, outdir:dir:details"
    envs = {
        "ignore_tcr": False,
        "tcell_selector": None,
        "indicator_genes": ["CD3E"],
        "kmeans": {"nstart": 25},
    }
    lang = config.lang.rscript
    script = "file://scripts/TCellSelection.R"
    plugin_opts = {"report": "file://reports/TCellSelection.svelte"}
    order = 5


class CloneHeterogeneity(Proc):
    """Clone heterogeneity in each cluster

    Envs:
        cases (type=json): Cases with keys as case names.
            Each case has arguments with keys:
            * cut: How to cut the clones by sizes
            * subsetting: How to subset the cells
            * design: Designed comparisons
    """
    input = "sobjfile:file"
    output = "outdir:dir:CloneHeterogeneity"
    lang = config.lang.rscript
    script = "file://scripts/CloneHeterogeneity.R"
    envs = {"cases": {}}
    plugin_opts = {
        "report": "file://reports/CloneHeterogeneity.svelte",
    }
    order = 11


class MetaMarkers(Proc):
    """Meta markers for different groups

    Envs:
        ncores (type=int): Number of cores to use
        cases (type=json): Cases with keys as case names and values.
            * filter: Focus on just subset of cells.
            * mutaters: Add new columns to the metadata.
            * group.by: Groups from which column to consider.
    """
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.meta_markers"
    lang = config.lang.rscript
    envs = {
        "ncores": 1,
        "cases": {},
    }
    script = "file://scripts/MetaMarkers.R"
    plugin_opts = {"report": "file://reports/MetaMarkers.svelte"}
    order = 10


class MarkersOverlapping(Proc):
    """Find the overlaping markers for multiple cases

    Envs:
        cases (type=json): Cases with keys as case names and values:
            * overlaps: The marker finder cases.
            * devpars: Devpars for the venn plot.

    """
    input = "mfdir:file"
    output = "outdir:dir:{{in.mfdir | stem}}.overlapping_markers"
    lang = config.lang.rscript
    envs = {
        "cases": {},
    }
    script = "file://scripts/MarkersOverlapping.R"
    plugin_opts = {"report": "file://reports/MarkersOverlapping.svelte"}
    order = 10
