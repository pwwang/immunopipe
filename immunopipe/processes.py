"""Process definition"""
from pipen.channel import expand_dir
from pipen_filters.filters import FILTERS
from biopipen.ns.misc import File2Proc
from biopipen.ns.tcr import (
    ImmunarchLoading,
    Immunarch,
    CloneResidency,
    Immunarch2VDJtools,
    VJUsage,
    TCRClustering,
    TCRClusteringStats,
)
from biopipen.ns.scrna import (
    SeuratPreparing,
    SeuratClustering,
    SeuratClusterStats,
    SeuratMetadataMutater as SeuratMetadataMutater_,
    MarkersFinder,
    CellsDistribution,
    ScFGSEA,
)
from biopipen.ns.scrna_metabolic import ScrnaMetabolic
from datar.all import tibble

from .args import config
from .inhouse import (
    SelectTCells,
    RadarPlots,
    CloneHeterogeneity,
    MetaMarkersForClones,
    MarkersOverlapping,
)

toml_dumps = FILTERS["toml_dumps"]
config = config or {}


class SampleInfo(File2Proc):
    """List sample information"""
    plugin_opts = {
        "report": "file://reports/SampleInfo.svelte",
        "report_toc": False,
    }


class ImmunarchLoading(ImmunarchLoading):
    """Load TCR data into immunarch object"""
    requires = SampleInfo


class SeuratPreparing(SeuratPreparing):
    requires = SampleInfo


class SeuratClusteringOfAllCells(SeuratClustering):
    """Cluster all cells"""
    requires = SeuratPreparing


class MarkersForClustersOfAllCells(MarkersFinder):
    """Markers for clusters of all cells"""
    requires = SeuratClusteringOfAllCells
    plugin_opts = {"report_order": 1}
    order = 4


class Immunarch(Immunarch):
    requires = ImmunarchLoading


class Immunarch2VDJtools(Immunarch2VDJtools):
    requires = ImmunarchLoading
    plugin_opts = {"args_hide": True}


class VJUsage(VJUsage):
    requires = Immunarch2VDJtools
    input_data = lambda ch: expand_dir(ch, pattern="*.txt")
    plugin_opts = {"report_toc": False}
    order = 2


# Start processes
STARTS = [SampleInfo]

if "SelectTCells" in config:
    class SelectTCells(SelectTCells):
        requires = [SeuratClusteringOfAllCells, ImmunarchLoading]

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


class SeuratClusterStats(SeuratClusterStats):
    requires = SeuratMetadataMutater
    order = 7


if "TCRClustering" in config or "TCRClusteringStats" in config:
    class TCRClustering(TCRClustering):
        requires = ImmunarchLoading

    class AttachTCRClusters2Seurat(SeuratMetadataMutater_):
        requires = SeuratClusteringOfTCells, TCRClustering
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1.rdsfile, metafile=ch2.clusterfile
        )

    class TCRClusteringStats(TCRClusteringStats):
        requires = TCRClustering


if "CellsDistribution" in config:
    class CellsDistribution(CellsDistribution):
        requires = (
            AttachTCRClusters2Seurat
            if "TCRClustering" in config or "TCRClusteringStats" in config
            else SeuratMetadataMutater
        )
        order = 8


if "CloneResidency" in config:
    class CloneResidency(CloneResidency):
        requires = ImmunarchLoading
        order = 3


if "CloneHeterogeneity" in config:
    class CloneHeterogeneity(CloneHeterogeneity):
        requires = SeuratMetadataMutater


if "RadarPlots" in config:
    class RadarPlots(RadarPlots):
        requires = SeuratMetadataMutater


if "ScFGSEA" in config:
    class ScFGSEA(ScFGSEA):
        requires = SeuratMetadataMutater


if "MarkersFinderForClones" in config:
    class MarkersFinderForClones(MarkersFinder):
        requires = SeuratMetadataMutater


if "MarkersOverlapping" in config:
    class MarkersOverlapping(MarkersOverlapping):
        requires = MarkersFinderForClones


if "MetaMarkersForClones" in config:
    class MetaMarkersForClones(MetaMarkersForClones):
        requires = SeuratMetadataMutater


if "METABOLIC" in config:
    MetabolicInputs = ScrnaMetabolic({"clustered": True}).starts[0]
    MetabolicInputs.order = 12
    MetabolicInputs.requires = SeuratMetadataMutater
    MetabolicInputs.input_data = lambda ch: tibble(
        metafile=ch.rdsfile,
        gmtfile=config["METABOLIC"]["gmtfile"],
        config=list(
            sorted(config["METABOLIC"]["cases"], key=lambda x: x.get("name"))
        ),
    )
