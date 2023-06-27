"""Process definition"""
from datar.tibble import tibble
from pipen.channel import expand_dir
from pipen.utils import mark
from pipen_annotate import annotate
from pipen_args import config
from pipen_board import from_pipen_board
from pipen_filters.filters import FILTERS

# biopipen processes
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
    CellTypeAnnotate as CellTypeAnnotate_,
    CellsDistribution,
    ScFGSEA,
)
from biopipen.ns.scrna_metabolic_landscape import ScrnaMetabolicLandscape

# inhouse processes
from .inhouse import (
    SelectTCells,
    RadarPlots,
    CloneHeterogeneity,
    MetaMarkersForClones,
    MarkersOverlapping,
)

toml_dumps = FILTERS["toml_dumps"]
from_board = from_pipen_board()


@annotate.format_doc(indent=1)
class SampleInfo(File2Proc):
    """Load and list sample information

    This process is just used to pass by input file and list the
    sample information in the report.

    Input:
        infile (required): {{Input.infile.help | indent: 12}}.
            The input file should have the following columns.
            * Sample: A unique id for each sample.
            * TCRDir: The directory for single-cell TCR data for this sample.
                Specifically, it should contain filtered_contig_annotations.csv
                or all_contig_annotations.csv from cellranger.
            * RNADir: The directory for single-cell RNA data for this sample.
                Specifically, it should be able to be read by
                `Seurat::Read10X()`.
                See also https://satijalab.org/seurat/reference/read10x.
            * Other columns are optional and will be treated as metadata for
                each sample.
    """
    plugin_opts = {
        "report": "file://reports/SampleInfo.svelte",
        "report_toc": False,
    }


class ImmunarchLoading(ImmunarchLoading):
    requires = SampleInfo


class SeuratPreparing(SeuratPreparing):
    requires = SampleInfo


@annotate.format_doc(indent=1)
class SeuratClusteringOfAllCells(SeuratClustering):
    """Cluster all cells using Seurat

    {{*Summary.long}}
    """
    requires = SeuratPreparing


@annotate.format_doc(indent=1)
class MarkersForClustersOfAllCells(MarkersFinder):
    """Find markers for clusters of all cells.

    If all your cells are T cells, the clustering will be performed on all
    T cells. `SeuratClusteringOfTCells` will be skipped.

    {{*Summary.long}}
    """
    requires = SeuratClusteringOfAllCells
    plugin_opts = {"report_order": 1}
    order = 4


class Immunarch(Immunarch):
    requires = ImmunarchLoading


@mark(board_config_hidden=True)
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

if "SelectTCells" in config or from_board:
    class SelectTCells(SelectTCells):
        requires = [SeuratClusteringOfAllCells, ImmunarchLoading]

    @annotate.format_doc(indent=2)
    class SeuratClusteringOfTCells(SeuratClustering):
        """Cluster the T cells selected by `SelectTCells`.

        If nothing is set for `SelectTCells` in the config file, meaning
        all cells are T cells, this process will be skipped.

        {{*Summary.long}}
        """
        requires = SelectTCells

    @annotate.format_doc(indent=2)
    class MarkersForClustersOfTCells(MarkersFinder):
        """Find markers for clusters of T cells

        This process will be skipped if all cells are T cells.

        {{*Summary.long}}
        """
        requires = SeuratClusteringOfTCells
        plugin_opts = {"report_order": 3}
        order = 6

else:
    SeuratClusteringOfTCells = SeuratClusteringOfAllCells


class CellTypeAnnotate(CellTypeAnnotate_):
    requires = SeuratClusteringOfTCells
    # Change the default to direct, which doesn't do any annotation
    envs = {"tool": "direct"}


@mark(board_config_hidden=True)
@annotate.format_doc(indent=1)
class SeuratMetadataMutater(SeuratMetadataMutater_):
    """Attach TCR clone information as meta columns to Seurat object

    {{*Summary.long}}
    """
    requires = CellTypeAnnotate, ImmunarchLoading
    input_data = lambda ch1, ch2: tibble(
        srtobj=ch1.outfile, metafile=ch2.metatxt
    )


class SeuratClusterStats(SeuratClusterStats):
    requires = SeuratMetadataMutater
    order = 7


if (
    "TCRClustering" in config
    or "TCRClusteringStats" in config
    or from_board
):
    @annotate.format_doc(indent=2)
    class TCRClustering(TCRClustering):
        """{{Summary.short}}

        You can disable this by remving the whole sections of
        TCRClustering and TCRClusteringStats in the config file.

        {{*Summary.long}}
        """
        requires = ImmunarchLoading

    @mark(board_config_hidden=True)
    @annotate.format_doc(indent=2)
    class AttachTCRClusters2Seurat(SeuratMetadataMutater_):
        """Attach TCR clusters as meta columns to Seurat object

        {{*Summary.long}}
        """
        requires = SeuratMetadataMutater, TCRClustering
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1.rdsfile, metafile=ch2.clusterfile
        )

    class TCRClusteringStats(TCRClusteringStats):
        requires = TCRClustering


if "CellsDistribution" in config or from_board:
    class CellsDistribution(CellsDistribution):
        requires = (
            AttachTCRClusters2Seurat
            if "TCRClustering" in config or "TCRClusteringStats" in config
            else SeuratMetadataMutater
        )
        order = 8


if "CloneResidency" in config or from_board:
    class CloneResidency(CloneResidency):
        requires = ImmunarchLoading
        order = 3


if "CloneHeterogeneity" in config or from_board:
    class CloneHeterogeneity(CloneHeterogeneity):
        requires = SeuratMetadataMutater


if "RadarPlots" in config or from_board:
    class RadarPlots(RadarPlots):
        requires = SeuratMetadataMutater


if "ScFGSEA" in config or from_board:
    class ScFGSEA(ScFGSEA):
        requires = SeuratMetadataMutater


if "MarkersFinderForClones" in config or from_board:
    class MarkersFinderForClones(MarkersFinder):
        requires = SeuratMetadataMutater


if "MarkersOverlapping" in config or from_board:
    class MarkersOverlapping(MarkersOverlapping):
        requires = MarkersFinderForClones


if "MetaMarkersForClones" in config or from_board:
    class MetaMarkersForClones(MetaMarkersForClones):
        requires = SeuratMetadataMutater


if "ScrnaMetabolicLandscape" in config or from_board:
    anno = annotate(ScrnaMetabolicLandscape)
    anno.Args.metafile.attrs["readonly"] = True
    anno.Args.metafile.attrs["hidden"] = True
    anno.Args.is_seurat.attrs["readonly"] = True
    anno.Args.is_seurat.attrs["hidden"] = True
    anno.Args.is_seurat.attrs["flag"] = True
    anno.Args.is_seurat.attrs["default"] = True
    anno.Args.is_seurat.attrs["value"] = True
    scrna_metabolic_landscape = ScrnaMetabolicLandscape(is_seurat=True)
    scrna_metabolic_landscape.p_input.order = 12
    scrna_metabolic_landscape.p_input.requires = SeuratMetadataMutater
