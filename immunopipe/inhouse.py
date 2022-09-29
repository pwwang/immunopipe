"""In house processes"""
from biopipen.core.proc import Proc
from biopipen.core.config import config


class SelectTCells(Proc):
    """Separate T and non-T cells"""
    input = "srtobj:file, immdata:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS, outdir:dir:details"
    envs = {"tcell_filter": "Clonotype_pct > .25", "indicator_gene": "CD3E"}
    lang = config.lang.rscript
    script = "file://scripts/SelectTCells.R"
    plugin_opts = {
        "report": "file://reports/SelectTCells.svelte",
        "report_toc": False,
        "report_order": 2,
    }
    order = 5


class RadarPlots(Proc):
    """Radar plots for cell proportion in different clusters"""
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.radar_plots"
    lang = config.lang.rscript
    script = "file://scripts/RadarPlots.R"
    envs = {"cases": {}}
    plugin_opts = {
        "report": "file://reports/RadarPlots.svelte",
        "report_toc": False,
    }


class CloneHeterogeneity(Proc):
    """Clone heterogeneity in each cluster"""
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


class MetaMarkersForClones(Proc):
    """Meta markers for different groups

    Envs:
        ncores: Number of cores to use
        cases: Cases with keys as case names and values -
            filter - Focus on just subset of cells
            mutaters - Add new columns to the metadata
            group.by - Groups from which column to consider

    """
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.meta_markers"
    lang = config.lang.rscript
    envs = {
        "ncores": 1,
        "cases": {},
    }
    script = "file://scripts/MetaMarkersForClones.R"
    plugin_opts = {"report": "file://reports/MetaMarkersForClones.svelte"}
    order = 10


class MarkersOverlapping(Proc):
    """Find the overlaping markers for multiple cases

    Envs:
        cases: Cases with keys as case names and values -
            overlaps - The marker finder cases
            devpars - Devpars for the venn plot

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
