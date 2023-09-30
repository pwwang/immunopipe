"""Process definition"""
from datar.tibble import tibble
from pipen.channel import expand_dir
from pipen.utils import mark
from pipen_annotate import annotate
from pipen_args import config
from pipen_board import from_pipen_board
from pipen_filters.filters import FILTERS

# biopipen processes
# from biopipen.ns.misc import File2Proc
from biopipen.ns.delim import SampleInfo as SampleInfo_
from biopipen.ns.tcr import (
    ImmunarchLoading as ImmunarchLoading_,
    Immunarch as Immunarch_,
    CloneResidency as CloneResidency_,
    Immunarch2VDJtools as Immunarch2VDJtools_,
    VJUsage as VJUsage_,
    TCRClustering as TCRClustering_,
    TCRClusteringStats as TCRClusteringStats_,
    CDR3AAPhyschem as CDR3AAPhyschem_,
    TESSA as TESSA_,
)
from biopipen.ns.scrna import (
    SeuratPreparing as SeuratPreparing_,
    SeuratClustering,
    SeuratClusterStats as SeuratClusterStats_,
    SeuratMetadataMutater as SeuratMetadataMutater_,
    MarkersFinder as MarkersFinder_,
    CellTypeAnnotation as CellTypeAnnotation_,
    CellsDistribution as CellsDistribution_,
    ScFGSEA as ScFGSEA_,
    TopExpressingGenes as TopExpressingGenes_,
    RadarPlots as RadarPlots_,
    ModuleScoreCalculator as ModuleScoreCalculator_,
    MetaMarkers as MetaMarkers_,
)
from biopipen.ns.scrna_metabolic_landscape import ScrnaMetabolicLandscape

# inhouse processes
from .inhouse import (
    TCellSelection as TCellSelection_,
    # CloneHeterogeneity,
    # MetaMarkers,
    # MarkersOverlapping,
)

toml_dumps = FILTERS["toml_dumps"]
from_board = from_pipen_board()


@annotate.format_doc(indent=1)
class SampleInfo(SampleInfo_):
    """{{Summary}}

    Input:
        infile (required): {{Input.infile.help | indent: 12}}.
            The input file should have the following columns.
            * Sample: A unique id for each sample.
            * TCRData: The directory for single-cell TCR data for this sample.
                Specifically, it should contain filtered_contig_annotations.csv
                or all_contig_annotations.csv from cellranger.
            * RNAData: The directory for single-cell RNA data for this sample.
                Specifically, it should be able to be read by
                `Seurat::Read10X()`.
                See also https://satijalab.org/seurat/reference/read10x.
            * Other columns are optional and will be treated as metadata for
                each sample.
    """
    envs = {"exclude_cols": "TCRData,RNAData"}


class ImmunarchLoading(ImmunarchLoading_):
    requires = SampleInfo


class Immunarch(Immunarch_):
    requires = ImmunarchLoading


@mark(board_config_hidden=True)
class Immunarch2VDJtools(Immunarch2VDJtools_):
    requires = ImmunarchLoading
    plugin_opts = {"args_hide": True}


class VJUsage(VJUsage_):
    requires = Immunarch2VDJtools
    input_data = lambda ch: expand_dir(ch, pattern="*.txt")
    plugin_opts = {"report_toc": False}
    order = 2


class SeuratPreparing(SeuratPreparing_):
    requires = SampleInfo


@annotate.format_doc(indent=1)
class SeuratClusteringOfAllCells(SeuratClustering):
    """Cluster all cells using Seurat

    {{*Summary.long}}
    """
    requires = SeuratPreparing


if "TCellSelection" in config or from_board:
    class TCellSelection(TCellSelection_):
        requires = [SeuratClusteringOfAllCells, ImmunarchLoading]

    @annotate.format_doc(indent=2)
    class SeuratClusteringOfTCells(SeuratClustering):
        """Cluster the T cells selected by `TCellSelection`.

        If nothing is set for `TCellSelection` in the config file, meaning
        all cells are T cells, this process will be skipped.

        {{*Summary.long}}
        """
        requires = TCellSelection
else:
    SeuratClusteringOfTCells = SeuratClusteringOfAllCells


class CellTypeAnnotation(CellTypeAnnotation_):
    requires = SeuratClusteringOfTCells
    # Change the default to direct, which doesn't do any annotation
    envs = {"tool": "direct"}


@annotate.format_doc(indent=1)
class ClusterMarkers(MarkersFinder_):
    """Markers for clusters of T cells.

    {{*Summary.long}}

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        ident-1 (hidden;readonly): {{Envs["ident-1"].help | indent: 12}}.
        ident-2 (hidden;readonly): {{Envs["ident-2"].help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
        prefix_each (hidden;readonly): {{Envs.prefix_each.help | indent: 12}}.
        section (hidden;readonly): {{Envs.section.help | indent: 12}}.
    """
    requires = CellTypeAnnotation
    envs = {"cases": {"Cluster": {}}}
    plugin_opts = {"report_order": 1}
    order = 4


@annotate.format_doc(indent=1)
class TopExpressingGenes(TopExpressingGenes_):
    """Top expressing genes for clusters of T cells.

    {{*Summary.long}}

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        ident (hidden;readonly): {{Envs.ident.help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
        prefix_each (hidden;readonly): {{Envs.prefix_each.help | indent: 12}}.
        section (hidden;readonly): {{Envs.section.help | indent: 12}}.
    """
    requires = CellTypeAnnotation
    envs = {"cases": {"Cluster": {}}}
    plugin_opts = {"report_order": 1}
    order = 4


if "ModuleScoreCalculator" in config or from_board:
    class ModuleScoreCalculator(ModuleScoreCalculator_):
        requires = CellTypeAnnotation

    CellTypeAnnotation = ModuleScoreCalculator


if "TESSA" in config or from_board:
    class TESSA(TESSA_):
        requires = ImmunarchLoading, CellTypeAnnotation
        input_data = lambda ch1, ch2: tibble(ch1.iloc[:, 0], ch2)

    CellTypeAnnotation = TESSA


# @mark(board_config_hidden=True)
@annotate.format_doc(indent=1)
class SeuratMetadataMutater(SeuratMetadataMutater_):
    """Attach TCR clone information as meta columns to Seurat object

    You may also use `envs.mutaters` to add new columns to the metadata.
    These columns can be used for downstream analysis.

    {{*Summary.long}}
    """
    requires = CellTypeAnnotation, ImmunarchLoading
    input_data = lambda ch1, ch2: tibble(
        srtobj=ch1.iloc[:, 0], metafile=ch2.metatxt
    )


if (
    "TCRClustering" in config
    or "TCRClusteringStats" in config
    or from_board
):
    @annotate.format_doc(indent=2)
    class TCRClustering(TCRClustering_):
        """{{Summary.short}}

        You can disable this by remving the whole sections of
        TCRClustering and TCRClusteringStats in the config file.

        {{*Summary.long}}
        """
        requires = ImmunarchLoading

    @mark(board_config_hidden=True)
    @annotate.format_doc(indent=2)
    class TCRClusters2Seurat(SeuratMetadataMutater_):
        """Attach TCR clusters as meta columns to Seurat object

        {{*Summary.long}}
        """
        requires = SeuratMetadataMutater, TCRClustering
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1.rdsfile, metafile=ch2.clusterfile
        )

    class TCRClusteringStats(TCRClusteringStats_):
        requires = TCRClustering

else:
    TCRClusters2Seurat = SeuratMetadataMutater


class SeuratClusterStats(SeuratClusterStats_):
    requires = TCRClusters2Seurat
    order = 7


if "CellsDistribution" in config or from_board:
    class CellsDistribution(CellsDistribution_):
        requires = TCRClusters2Seurat
        order = 8


if "CloneResidency" in config or from_board:
    class CloneResidency(CloneResidency_):
        requires = ImmunarchLoading
        order = 3


# if "CloneHeterogeneity" in config or from_board:
#     class CloneHeterogeneity(CloneHeterogeneity):
#         requires = SeuratMetadataMutater


if "RadarPlots" in config or from_board:
    class RadarPlots(RadarPlots_):
        requires = TCRClusters2Seurat


if "ScFGSEA" in config or from_board:
    class ScFGSEA(ScFGSEA_):
        requires = TCRClusters2Seurat


if "MarkersFinder" in config or from_board:
    class MarkersFinder(MarkersFinder_):
        requires = TCRClusters2Seurat


# if "MarkersOverlapping" in config or from_board:
#     class MarkersOverlapping(MarkersOverlapping):
#         requires = MarkersFinder


if "MetaMarkers" in config or from_board:
    class MetaMarkers(MetaMarkers_):
        requires = TCRClusters2Seurat


if "CDR3AAPhyschem" in config or from_board:
    class CDR3AAPhyschem(CDR3AAPhyschem_):
        requires = ImmunarchLoading, TCRClusters2Seurat
        input_data = lambda ch1, ch2: tibble(
            immdata=ch1.rdsfile,
            srtobj=ch2.rdsfile,
        )


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
    scrna_metabolic_landscape.p_input.requires = TCRClusters2Seurat
    scrna_metabolic_landscape.p_input.order = 99
