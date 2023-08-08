"""In house processes"""
from biopipen.core.proc import Proc
from biopipen.core.config import config


class TCellSelection(Proc):
    """Separate T and non-T cells and select T cells

    If all of your cells are T cells, you can skip this process by removing
    the whole section `[TCellSelection]` from the configuration file.

    In such a case, `SeuratClusteringOfAllCells` will be clustering all
    T cells and `SeuratClusteringOfTCells` is skipped.

    Input:
        srtobj: Seurat object file in RDS
        immdata: Immunarch data file in RDS

    Output:
        rdsfile: Seurat object file in RDS
        outdir: Output directory with details

    Envs:
        tcell_indicator: Use `Clonotype_Pct` directly to determine T cells
            This will be passed to `immdata |> mutate(is_TCell = ...)`
            For example, `Clonotype_Pct > 0.25`.
        indicator_genes (list): A list of indicator genes whose expression
            values and clonotype percentage will be used to determine T cells.
            A kmeans clustering will be performed on those values with K=2.
    """
    input = "srtobj:file, immdata:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS, outdir:dir:details"
    envs = {
        "tcell_indicator": None,
        "indicator_genes": ["CD3E"],
    }
    lang = config.lang.rscript
    script = "file://scripts/TCellSelection.R"
    plugin_opts = {
        "report": "file://reports/TCellSelection.svelte",
        # "report_toc": False,
        "report_order": 2,
    }
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
        "report_order": 20,
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
