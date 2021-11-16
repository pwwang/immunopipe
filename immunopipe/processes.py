"""Process definition"""
from itertools import islice
from pipen import Proc
from pipen.channel import expand_dir
from pipen_filters.filters import FILTERS
from biopipen.core.filters import filtermanager
from biopipen.namespaces.misc import File2Proc, Config2File
from biopipen.namespaces.csv import BindRows
from biopipen.namespaces.tcr import (
    ImmunarchLoading,
    ImmunarchBasic,
    ImmunarchAdvanced,
    ImmunarchFilter,
    CloneResidency,
    Immunarch2VDJtools,
    VJUsage,
)
from biopipen.namespaces.scrna import (
    SeuratPreparing,
    SeuratClustering,
    MarkersFinder,
    DimPlots,
    GeneExpressionInvistigation,
)
from datar.all import (
    f,
    bind_cols,
    select,
    rename,
    unite,
    tibble,
    everything,
    rowwise,
    mutate,
    strsplit,
    flatten,
)

from .args import config

starts = []


def chunk_list(array, lens):
    start = 0
    out = []
    for length in lens:
        out.append(array[start : (start + length)])
        start += length
    return out


class SampleInfo(File2Proc):
    """List sample information"""

    plugin_opts = {
        "report": "file://reports/SampleInfo.svelte",
        "report_toc": False,
    }


GeneList = Proc.from_proc(File2Proc)
starts.append(SampleInfo)
starts.append(GeneList)

ImmunarchLoading = Proc.from_proc(
    ImmunarchLoading, requires=SampleInfo, plugin_opts={"args_hide": True}
)

SeuratPreparing = Proc.from_proc(
    SeuratPreparing,
    requires=SampleInfo,
    plugin_opts={"args_hide": True},
)

SeuratClusteringOfAllCells = Proc.from_proc(
    SeuratClustering,
    requires=SeuratPreparing,
    plugin_opts={"args_hide": True},
)

MarkersForClustersOfAllCells = Proc.from_proc(
    MarkersFinder,
    requires=SeuratClusteringOfAllCells,
    input_data=lambda ch: tibble(
        ch,
        "Markers for clusters of all cells",
    ),
    envs={"cases": "ident"},
    plugin_opts={
        "report_order": 1,
        "args_hide": True,
    },
)


class SelectTCells(Proc):
    """Separate T and non-T cells"""

    input = "srtobj:file, immdata:file"
    requires = [SeuratClusteringOfAllCells, ImmunarchLoading]
    input_data = lambda ch1, ch2: tibble(
        select(ch1, 1),
        ch2,
        _name_repair="minimal",
    )
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS, outdir:dir:details"
    envs = {"tcell_filter": "Clonotype_pct > .25"}
    lang = MarkersForClustersOfAllCells.lang
    script = "file://scripts/SelectTCells.R"
    plugin_opts = {
        "report": "file://reports/SelectTCells.svelte",
        "report_toc": False,
        "report_order": 2,
        "args_hide": True,
    }


SeuratClusteringOfTCells = Proc.from_proc(
    SeuratClustering,
    requires=SelectTCells,
    plugin_opts={"args_hide": True},
)

MarkersForClustersOfTCells = Proc.from_proc(
    MarkersFinder,
    requires=SeuratClusteringOfTCells,
    input_data=lambda ch: tibble(ch, "Markers for clusters of T cells"),
    envs={"cases": "ident"},
    plugin_opts={"report_order": 3, "args_hide": True},
)

ImmunarchBasic = Proc.from_proc(
    ImmunarchBasic,
    requires=ImmunarchLoading,
    plugin_opts={"args_hide": True},
)

ImmunarchAdvanced = Proc.from_proc(
    ImmunarchAdvanced,
    requires=ImmunarchLoading,
    plugin_opts={"args_hide": True},
)

CloneResidency = Proc.from_proc(
    CloneResidency,
    requires=ImmunarchLoading,
    plugin_opts={"args_hide": True},
)

Immunarch2VDJtools = Proc.from_proc(
    Immunarch2VDJtools,
    requires=ImmunarchLoading,
    plugin_opts={"args_hide": True},
)

VJUsage = Proc.from_proc(
    VJUsage,
    requires=Immunarch2VDJtools,
    input_data=lambda ch: expand_dir(ch, pattern="*.txt"),
    plugin_opts={"report_toc": False, "args_hide": True},
)

if "RADAR_PLOTS" in config:
    RadarPlotsConfig = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](conf["filters"][key])
            for conf in config.RADAR_PLOTS
            for key in conf["filters"]
        ],
    )
    starts.append(RadarPlotsConfig)

    RadarPlotsFilter = Proc.from_proc(
        ImmunarchFilter,
        requires=[ImmunarchLoading, RadarPlotsConfig],
        input_data=lambda ch1, ch2: tibble(
            ch1,
            ch2,
            [key for conf in config.RADAR_PLOTS for key in conf["filters"]],
        ),
    )

    class RadarPlots(Proc):
        """Radar plots for cell proportion in different clusters"""

        input = [
            "srtobj:file",
            "groupfiles:files",
            "name:var",
            "direction:var",
            "breaks:var",
        ]
        requires = [SeuratClusteringOfTCells, RadarPlotsFilter]
        input_data = lambda ch1, ch2: tibble(
            select(ch1, 1),
            chunk_list(
                flatten(select(ch2, 2)),
                [len(conf["filters"]) for conf in config.RADAR_PLOTS],
            ),
            [conf.name for conf in config.RADAR_PLOTS],
            [
                conf.get("direction", "intra-cluster")
                for conf in config.RADAR_PLOTS
            ],
            [
                ",".join(str(brk) for brk in conf.get("breaks", [0, 50, 100]))
                for conf in config.RADAR_PLOTS
            ],
        )
        output = "outfile:file:{{in.name | slugify}}.radar.png"
        lang = SeuratClusteringOfTCells.lang
        script = "file://scripts/RadarPlots.R"
        plugin_opts = {
            "report": "file://reports/RadarPlots.svelte",
            "report_toc": False,
        }
        template_opts = {"filters": filtermanager.filters.copy()}


if "MARKERS_FINDER" in config:
    markers_groups = list(
        [case for case in config.MARKERS_FINDER.filters[0] if case != "name"]
    )
    MarkersFinderClonesFilterConfig = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](conf[markers_groups[0]])
            for conf in config.MARKERS_FINDER.filters
        ],
    )
    MarkersFinderRestClonesFilterConfig = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](conf[markers_groups[1]])
            for conf in config.MARKERS_FINDER.filters
        ],
    )
    starts.append(MarkersFinderClonesFilterConfig)
    starts.append(MarkersFinderRestClonesFilterConfig)

    MarkersFinderClonesFilter = Proc.from_proc(
        ImmunarchFilter,
        requires=[ImmunarchLoading, MarkersFinderClonesFilterConfig],
        input_data=lambda ch1, ch2: tibble(ch1, ch2, name=markers_groups[0]),
    )

    MarkersFinderRestClonesFilter = Proc.from_proc(
        ImmunarchFilter,
        requires=[ImmunarchLoading, MarkersFinderRestClonesFilterConfig],
        input_data=lambda ch1, ch2: tibble(ch1, ch2, name=markers_groups[1]),
    )

    MarkersFinderBindGroups = Proc.from_proc(
        BindRows,
        requires=[MarkersFinderClonesFilter, MarkersFinderRestClonesFilter],
        input_data=lambda ch1, ch2: (
            ch1
            >> rename(groupfile0=f.groupfile, o=f.outfile)  # suppress warning
            >> bind_cols(ch2)
            >> select(f.groupfile0, f.groupfile)
            >> unite(f.infiles, everything(), sep=" ||| ")
            >> rowwise()
            >> mutate(infiles=strsplit(f.infiles, " ||| ", fixed=True))
        ),
    )
    MarkersFinder = Proc.from_proc(
        MarkersFinder,
        requires=[SeuratPreparing, MarkersFinderBindGroups],
        input_data=lambda ch1, ch2: tibble(
            ch1,
            ch2,
            [filt["name"] for filt in config.MARKERS_FINDER.filters],
        ),
        envs={"cases": {markers_groups[0]: {"ident.2": markers_groups[1]}}},
    )

if "GENE_EXPR_INVESTIGATION_CLUSTERS" in config:
    # Do the filtering
    GeneExprInvestigationClustersConfig = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](conf)
            for conf in config.GENE_EXPR_INVESTIGATION_CLUSTERS
        ],
    )
    starts.append(GeneExprInvestigationClustersConfig)
    from datar.all import t

    GeneExpressionInvistigationClusters = Proc.from_proc(
        GeneExpressionInvistigation,
        requires=[
            SeuratClusteringOfTCells,
            GeneList,
            GeneExprInvestigationClustersConfig,
        ],
        input_data=lambda ch1, ch2, ch3: tibble(ch1, t(ch2), ch3),
    )

if "DIM_PLOTS" in config:
    DimPlotsConfig = Proc.from_proc(
        Config2File,
        input_data=[FILTERS["toml_dumps"](conf) for conf in config.DIM_PLOTS],
    )
    starts.append(DimPlotsConfig)
    DimPlots = Proc.from_proc(
        DimPlots,
        requires=[SeuratClusteringOfTCells, DimPlotsConfig],
    )
