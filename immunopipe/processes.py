"""Process definition"""
import numpy
from pipen import Proc
from pipen.channel import expand_dir
from pipen_filters.filters import FILTERS
from biopipen.core.filters import filtermanager
from biopipen.namespaces.misc import File2Proc, Config2File
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
    GeneExpressionInvestigation,
)
from biopipen.namespaces.scrna_metabolic import (
    build_processes,
)
from datar.all import (
    select,
    tibble,
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
    cache="force",
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
    cache="force",
    requires=SeuratClusteringOfAllCells,
    input_data=lambda ch: tibble(
        ch,
        None,
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
    cache="force",
    requires=SeuratClusteringOfTCells,
    input_data=lambda ch: tibble(ch, None, "Markers for clusters of T cells"),
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
            FILTERS["toml_dumps"](conf.filters) for conf in config.RADAR_PLOTS
        ],
    )
    starts.append(RadarPlotsConfig)

    RadarPlotsFilter = Proc.from_proc(
        ImmunarchFilter,
        requires=[ImmunarchLoading, RadarPlotsConfig],
        envs={"merge": True},
    )

    class RadarPlots(Proc):
        """Radar plots for cell proportion in different clusters"""

        input = [
            "srtobj:file",
            "groupfile:file",
            "name:var",
            "direction:var",
            "breaks:var",
        ]
        requires = [SeuratClusteringOfTCells, RadarPlotsFilter]
        input_data = lambda ch1, ch2: tibble(
            select(ch1, 1),
            select(ch2, 2),
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
    ImmunarchFilterConfig = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](mf_config.subsetting)
            for mf_config in config.MARKERS_FINDER
        ],
    )
    MarkersFinderCases = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](
                {
                    name: {"ident.1": design[0], "ident.2": design[1]}
                    for name, design in sorted(mf_config.design.items())
                }
            )
            for mf_config in config.MARKERS_FINDER
        ],
    )
    MetaMarkersConfig = Proc.from_proc(
        Config2File,
        input_data=[
            FILTERS["toml_dumps"](
                {
                    "name": mf_config.name,
                    "meta": mf_config.meta,
                    "overlap": mf_config.overlap,
                }
            )
            for mf_config in config.MARKERS_FINDER
        ],
    )
    starts.append(ImmunarchFilterConfig)
    starts.append(MarkersFinderCases)
    starts.append(MetaMarkersConfig)

    MarkersFinderClonesFilter = Proc.from_proc(
        ImmunarchFilter,
        requires=[ImmunarchLoading, ImmunarchFilterConfig],
        envs={"merge": True},
    )

    MarkersFinderClones = Proc.from_proc(
        MarkersFinder,
        cache="force",
        requires=[
            SeuratPreparing,
            MarkersFinderClonesFilter,
            MarkersFinderCases,
        ],
        input_data=lambda ch1, ch2, ch3: tibble(
            srtobj=ch1,
            groupfile=select(ch2, 2),
            casefile=ch3,
            name=[mf_config.name for mf_config in config.MARKERS_FINDER],
        ),
    )

    class MetaMarkersForClones(Proc):
        """Meta markers for different groups"""
        cache = "force"
        requires = (
            SeuratPreparing,
            MarkersFinderClonesFilter,
            MarkersFinderClones,
            MetaMarkersConfig,
        )
        input = "srtobj:file, groupfile:file, markersdir:file, configfile:file"
        input_data = lambda ch1, ch2, ch3, ch4: tibble(
            srtobj=ch1,
            groupfile=ch2.iloc[:, 1],
            markersdir=ch3,
            configfile=ch4,
        )
        output = "outdir:dir:MetaMarkers"
        lang = SeuratPreparing.lang
        envs = {
            "ncores": 1,
            "min_cells": 5,
            "venn_devpars": {"res": 100, "height": 1000, "width": 1200},
            "upset_devpars": {"res": 100, "height": 1000, "width": 1000},
        }
        script = "file://scripts/MetaMarkersForClones.R"
        plugin_opts = {
            "report": "file://reports/MetaMarkersForClones.svelte"
        }
        template_opts = SeuratPreparing.template_opts


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

    GeneExpressionInvestigationClusters = Proc.from_proc(
        GeneExpressionInvestigation,
        requires=[
            SeuratClusteringOfTCells,
            GeneList,
            GeneExprInvestigationClustersConfig,
        ],
        input_data=lambda ch1, ch2, ch3: tibble(ch1, [flatten(ch2)], ch3),
        envs={"gopts": {"header": False, "sep": "\t", "row.names": None}},
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

if "METABOLIC" in config:
    if any(
        conf_case.get("subset_using") == "immunarch"
        for conf_case in sorted(
            config["METABOLIC"]["cases"], key=lambda x: x["name"]
        )
    ):
        MetabolicSubsetFilterConfig = Proc.from_proc(
            Config2File,
            input_data=[
                FILTERS["toml_dumps"](conf_case["subsetting"])
                for conf_case in sorted(
                    config["METABOLIC"]["cases"], key=lambda x: x["name"]
                )
                if conf_case.get("subset_using") == "immunarch"
            ],
        )
        MetabolicSubsetFilter = Proc.from_proc(
            ImmunarchFilter,
            requires=[ImmunarchLoading, MetabolicSubsetFilterConfig],
            envs={"merge": True},
        )
        starts.append(MetabolicSubsetFilterConfig)
    else:
        MetabolicSubsetFilter = None

    MetabolicInputs = build_processes({"clustered": True})
    if MetabolicSubsetFilter:

        def metabolic_input(ch1, ch2):
            subsetfile = numpy.array(
                [None] * len(config["METABOLIC"]["cases"])
            )
            subsetfile[
                [
                    conf_case.get("subset_using") == "immunarch"
                    for conf_case in sorted(
                        config["METABOLIC"]["cases"], key=lambda x: x["name"]
                    )
                ]
            ] = ch2.groupfile
            return tibble(
                metafile=ch1.rdsfile,
                subsetfile=subsetfile,
                groupfile=None,
                gmtfile=config["METABOLIC"]["gmtfile"],
                config=[
                    FILTERS["toml_dumps"](conf_case)
                    for conf_case in sorted(
                        config["METABOLIC"]["cases"], key=lambda x: x["name"]
                    )
                ],
            )

        MetabolicInputs.requires = (
            SeuratClusteringOfTCells,
            MetabolicSubsetFilter,
        )
        MetabolicInputs.input_data = metabolic_input
    else:
        MetabolicInputs.requires = SeuratClusteringOfTCells
        MetabolicInputs.input_data = lambda ch: tibble(
            metafile=ch.rdsfile,
            subsetfile=None,
            groupfile=None,
            gmtfile=config["METABOLIC"]["gmtfile"],
            config=[
                FILTERS["toml_dumps"](conf_case)
                for conf_case in sorted(
                    config["METABOLIC"]["cases"], key=lambda x: x["name"]
                )
            ],
        )
