"""Process definition"""
import numpy
from pipen import Proc
from pipen.channel import expand_dir
from pipen_filters.filters import FILTERS
from biopipen.core.filters import filtermanager
from biopipen.namespaces.misc import File2Proc, Config2File
from biopipen.namespaces.tcr import (
    ImmunarchLoading as ImmunarchLoading_,
    Immunarch as Immunarch_,
    ImmunarchFilter,
    CloneResidency as CloneResidency_,
    Immunarch2VDJtools as Immunarch2VDJtools_,
    VJUsage as VJUsage_,
)
from biopipen.namespaces.scrna import (
    SeuratPreparing as SeuratPreparing_,
    SeuratClustering,
    SeuratMetadataMutater as SeuratMetadataMutater_,
    MarkersFinder,
    DimPlots as DimPlots_,
    GeneExpressionInvestigation,
)
from biopipen.namespaces.scrna_metabolic import build_processes
from datar.all import f, select, tibble

from .args import config
from .utils import chunk_list


class SampleInfo(File2Proc):
    """List sample information"""

    plugin_opts = {
        "report": "file://reports/SampleInfo.svelte",
        "report_toc": False,
    }


class ImmunarchLoading(ImmunarchLoading_):
    """Load TCR data into immunarch object"""
    requires = SampleInfo


class SeuratPreparing(SeuratPreparing_):
    requires = SampleInfo


class SeuratClusteringOfAllCells(SeuratClustering):
    """Cluster all cells"""
    requires = SeuratPreparing


class MarkersForClustersOfAllCells(MarkersFinder):
    """Markers for clusters of all cells"""
    requires = SeuratClusteringOfAllCells
    plugin_opts = {"report_order": 1}
    order = 4


class Immunarch(Immunarch_):
    requires = ImmunarchLoading


class Immunarch2VDJtools(Immunarch2VDJtools_):
    requires = ImmunarchLoading
    plugin_opts = {"args_hide": True}


class VJUsage(VJUsage_):
    requires = Immunarch2VDJtools
    input_data = lambda ch: expand_dir(ch, pattern="*.txt")
    plugin_opts = {"report_toc": False}
    order = 2


# Start processes
STARTS = [SampleInfo]

if config.get("SelectTCells", True):

    class SelectTCells(Proc):
        """Separate T and non-T cells"""
        input = "srtobj:file, immdata:file"
        requires = [SeuratClusteringOfAllCells, ImmunarchLoading]
        output = "rdsfile:file:{{in.srtobj | stem}}.RDS, outdir:dir:details"
        envs = {"tcell_filter": "Clonotype_pct > .25", "indicator_gene": "CD3E"}
        lang = MarkersForClustersOfAllCells.lang
        script = "file://scripts/SelectTCells.R"
        plugin_opts = {
            "report": "file://reports/SelectTCells.svelte",
            "report_toc": False,
            "report_order": 2,
        }
        order = 5

    class SeuratClusteringOfTCells(SeuratClustering):
        requires = SelectTCells

    class MarkersForClustersOfTCells(MarkersFinder):
        requires = SeuratClusteringOfTCells
        plugin_opts = {"report_order": 3}
        order = 6

else:
    SeuratClusteringOfTCells = SeuratClusteringOfAllCells


class SeuratMetadataMutater(SeuratMetadataMutater_):
    requires = SeuratClusteringOfTCells, ImmunarchLoading
    input_data = lambda ch1, ch2: tibble(
        srtobj=ch1.rdsfile, metafile=ch2.metatxt
    )


class DimPlots(DimPlots_):
    requires = SeuratMetadataMutater
    order = 7


if "CloneResidency" in config:
    class CloneResidency(CloneResidency_):
        requires = ImmunarchLoading
        order = 3


if "CloneHeterogeneity" in config:

    class CloneHeterogeneity(Proc):
        """Clone heterogeneity in each cluster"""
        input = "sobjfile:file"
        requires = SeuratMetadataMutater
        output = "outdir:dir:CloneHeterogeneity"
        lang = DimPlots.lang
        script = "file://scripts/CloneHeterogeneity.R"
        envs = {"cases": {}}
        template_opts = {"filters": filtermanager.filters.copy()}
        plugin_opts = {
            "report": "file://reports/CloneHeterogeneity.svelte",
            "report_order": 20,
        }
        order = 11


if "RADAR_PLOTS" in config:

    class RadarPlotsConfig(Config2File):
        input_data = [
            FILTERS["toml_dumps"](conf) for conf in config.RADAR_PLOTS
        ]

    class RadarPlotsFilter(ImmunarchFilter):
        requires = [ImmunarchLoading, RadarPlotsConfig]
        input_data = lambda ch1, ch2: tibble(
            immdata=select(ch1, 0),
            filterfile=ch2,
        )

    class RadarPlots(Proc):
        """Radar plots for cell proportion in different clusters"""

        input = [
            "srtobj:file",
            "groupfile:file",
            "direction:var",
            "breaks:var",
        ]
        requires = [SeuratClusteringOfTCells, RadarPlotsFilter]
        input_data = lambda ch1, ch2: tibble(
            select(ch1, 0),
            select(ch2, 1),
            [
                conf.get("direction", "intra-cluster")
                for conf in config.RADAR_PLOTS
            ],
            [
                ",".join(str(brk) for brk in conf.get("breaks", [0, 50, 100]))
                for conf in config.RADAR_PLOTS
            ],
            _name_repair="minimal",
        )
        output = "outfile:file:{{in.groupfile | stem | slugify}}.radar.png"
        lang = SeuratClusteringOfTCells.lang
        script = "file://scripts/RadarPlots.R"
        plugin_opts = {
            "report": "file://reports/RadarPlots.svelte",
            "report_toc": False,
        }
        template_opts = {"filters": filtermanager.filters.copy()}

    STARTS.append(RadarPlotsConfig)


if "MARKERS_FINDER" in config:

    class ImmunarchFilterConfig(Config2File):
        input_data = [
            FILTERS["toml_dumps"](mf_config)
            for mf_config in config.MARKERS_FINDER
        ]

    class MarkersFinderClonesFilter(ImmunarchFilter):
        requires = [ImmunarchLoading, ImmunarchFilterConfig]
        input_data = lambda ch1, ch2: tibble(
            immdata=select(ch1, 0),
            filterfile=ch2,
        )
        order = 9

    class ApplyFiltersToSeurat(SeuratMetadataMutater_):
        requires = SeuratPreparing, MarkersFinderClonesFilter
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1,
            metafile=select(ch2, 1)
        )

    class MarkersFinderCases(Config2File):
        input_data = [
            FILTERS["toml_dumps"](
                dict(
                    name=mf_config.name,
                    cases={
                        name: {
                            "ident.1": design[0],
                            "ident.2": design[1],
                            "group.by": mf_config.name,
                        }
                        for name, design in sorted(mf_config.design.items())
                    }
                )
            )
            for mf_config in config.MARKERS_FINDER
        ]

    class MarkersFinderClones(MarkersFinder):
        requires = [
            ApplyFiltersToSeurat,
            MarkersFinderCases,
        ]
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1,
            casefile=ch2,
        )


    class MetaMarkersConfig(Config2File):
        input_data = [
            FILTERS["toml_dumps"](
                {
                    "name": mf_config.name,
                    "meta": mf_config.get("meta"),
                    "overlap": mf_config.get("overlap"),
                }
            )
            for mf_config in config.MARKERS_FINDER
            if "meta" in mf_config or "overlap" in mf_config
        ]

    class MetaMarkersForClones(Proc):
        """Meta markers for different groups"""

        # cache = "force"
        if MetaMarkersConfig.input_data:
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
        plugin_opts = {"report": "file://reports/MetaMarkersForClones.svelte"}
        template_opts = SeuratPreparing.template_opts
        order = 10

    STARTS.append(ImmunarchFilterConfig)
    STARTS.append(MarkersFinderCases)
    if MetaMarkersConfig.input_data:
        STARTS.append(MetaMarkersConfig)


if "GENE_EXPR_INVESTIGATION_CLUSTERS" in config:
    # Do the filtering
    class GeneExprInvestigationClustersConfig(Config2File):
        input_data = [
            FILTERS["toml_dumps"](conf)
            for conf in config.GENE_EXPR_INVESTIGATION_CLUSTERS
        ]

    class GeneExpressionInvestigationClusters(GeneExpressionInvestigation):
        requires = [
            SeuratClusteringOfTCells,
            GeneExprInvestigationClustersConfig,
        ]
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1,
            genefile=[
                conf.genefile
                for conf in config.GENE_EXPR_INVESTIGATION_CLUSTERS
            ],
            configfile=ch2,
        )
        order = 8

    STARTS.append(GeneExprInvestigationClustersConfig)


if "METABOLIC" in config:

    MetabolicInputs = build_processes({"clustered": True})
    MetabolicInputs.order = 12
    MetabolicInputs.requires = SeuratClusteringOfTCells
    MetabolicInputs.input_data = lambda ch: tibble(
        metafile=ch.rdsfile,
        gmtfile=config["METABOLIC"]["gmtfile"],
        config=list(
            sorted(config["METABOLIC"]["cases"], key=lambda x: x.get("name"))
        ),
    )
