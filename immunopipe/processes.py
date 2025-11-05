"""Process definition"""

from __future__ import annotations
from typing import Type, Sequence, Callable

from pipen.utils import is_loading_pipeline
from pipen_annotate import annotate
from pipen_filters.filters import FILTERS

# biopipen processes
from biopipen.core.config import config as biopipen_config
from biopipen.core.proc import Proc
from biopipen.ns.delim import SampleInfo as SampleInfo_
from biopipen.ns.tcr import (
    ScRepLoading as ScRepLoading_,
    # Immunarch as Immunarch_,
    # CloneResidency as CloneResidency_,
    TCRClustering as TCRClustering_,
    # TCRClusterStats as TCRClusterStats_,
    CDR3AAPhyschem as CDR3AAPhyschem_,
    TESSA as TESSA_,
    ScRepCombiningExpression as ScRepCombiningExpression_,
    ClonalStats as ClonalStats_,
)
from biopipen.ns.scrna import (
    SeuratPreparing as SeuratPreparing_,
    SeuratClustering as SeuratClustering_,
    SeuratSubClustering as SeuratSubClustering_,
    SeuratMap2Ref as SeuratMap2Ref_,
    SeuratClusterStats as SeuratClusterStats_,
    # SeuratMetadataMutater as SeuratMetadataMutater_,
    MarkersFinder as MarkersFinder_,
    CellTypeAnnotation as CellTypeAnnotation_,
    # CellsDistribution as CellsDistribution_,
    ScFGSEA as ScFGSEA_,
    TopExpressingGenes as TopExpressingGenes_,
    ModuleScoreCalculator as ModuleScoreCalculator_,
    CellCellCommunication as CellCellCommunication_,
    CellCellCommunicationPlots as CellCellCommunicationPlots_,
    PseudoBulkDEG as PseudoBulkDEG_,
)
from biopipen.ns.scrna_metabolic_landscape import ScrnaMetabolicLandscape

# inhouse processes
from .inhouse import (
    # TCellSelection as TCellSelection_,
    TOrBCellSelection as TOrBCellSelection_,
)
from .validate_config import validate_config

toml_dumps = FILTERS["toml_dumps"]
just_loading = is_loading_pipeline("help", "-h", "--help", "-h+", "--help+")
config = validate_config()

# https://pwwang.github.io/immunopipe/latest/
TEST_OUTPUT_BASEURL = "https://raw.githubusercontent.com/pwwang/immunopipe/tests-output"
start_processes = []


def when(
    condition: bool,
    requires: Type[Proc] | Sequence[Type[Proc]] | None = None,
) -> Callable[[Type[Proc]], Type[Proc] | None]:
    """Decorator to conditionally define the processes

    Args:
        condition: The condition to check
            If False, None is returned and the process is not defined.
        requires: The requirements of the process

    Returns:
        The decorator function
    """

    def decorator(cls: Type[Proc]) -> Type[Proc]:
        if condition or just_loading:
            if requires is not None:
                cls.requires = requires
            else:
                start_processes.append(cls)
            return cls

        return None

    return decorator


# Either has both RNA and VDJ data, or just VDJ data when LoadingRNAFromSeurat is used
@when("SampleInfo" in config or "LoadingRNAFromSeurat" not in config)
@annotate.format_doc(vars={"output_baseurl": TEST_OUTPUT_BASEURL})
class SampleInfo(SampleInfo_):
    __doc__ = """{{Summary}}

    This process is the entrance of the pipeline. It just pass by input file and list
    the sample information in the report.

    To specify the input file in the configuration file, use the following

    ```toml
    [SampleInfo.in]
    infile = [ "path/to/sample_info.txt" ]
    ```

    Or with `pipen-board`, find the `SampleInfo` process and click the `Edit` button.
    Then you can specify the input file here

    ![infile](images/SampleInfo-infile.png)

    Multiple input files are supported by the underlying pipeline framework. However,
    we recommend to run it with a different pipeline instance with configuration files.

    For the content of the input file, please see details
    [here](../preparing-input.md#metadata).

    You can add some columns to the input file while doing the statistics or you can
    even pass them on to the next processes. See `envs.mutaters` and
    `envs.save_mutated`.
    But if you are adding a factor (categorical) column with desired levels, the order
    can't be guaranteed, because we are saving them to a text file, where we can't
    guarantee the order of the levels. If you want to add a factor column with desired
    levels, you can set `envs.mutaters` of the `SeuratPreparing` process to mutate the
    column.

    Once the pipeline is finished, you can see the sample information in the report

    ![report](images/SampleInfo-report.png)

    Note that the required `RNAData` (if not loaded from a Seurat object) and
    `TCRData`/`BCRData` columns are not shown in the report.
    They are used to specify the paths of the `scRNA-seq` and `scTCR-seq`/`scBCR-seq`
    data, respectively.
    Also note that when `RNAData` is loaded from a Seurat object (specified in the
    `LoadingRNAFromSeurat` process), the metadata provided in this process will not be
    integrated into the Seurat object in the downstream processes. To incoporate
    these meta information into the Seurat object, please provide them in the
    Seurat object itself or use the `envs.mutaters` of the `SeuratPreparing` process
    to mutate the metadata of the Seurat object. But the meta information provided in
    this process can still be used in the statistics and plots in the report.

    You may also perform some statistics on the sample information, for example,
    number of samples per group. See next section for details.

    /// Tip
    This is the start process of the pipeline. Once you change the parameters for
    this process, the whole pipeline will be re-run.

    If you just want to change the parameters for the statistics, and use the
    cached (previous) results for other processes, you can set `cache` at
    pipeline level to `"force"` to force the pipeline to use the cached results
    and `cache` of `SampleInfo` to `false` to force the pipeline to re-run the
    `SampleInfo` process only.

    ```toml
    cache = "force"

    [SampleInfo]
    cache = false
    ```
    ///

    Examples:
        ### Example data

        | Sample | Age | Sex | Diagnosis |
        |--------|-----|-----|-----------|
        | C1     | 62  | F   | Colitis   |
        | C2     | 71.2| F   | Colitis   |
        | C3     | 56.2| M   | Colitis   |
        | C4     | 61.5| M   | Colitis   |
        | C5     | 72.8| M   | Colitis   |
        | C6     | 78.4| M   | Colitis   |
        | C7     | 61.6| F   | Colitis   |
        | C8     | 49.5| F   | Colitis   |
        | NC1    | 43.6| M   | NoColitis |
        | NC2    | 68.1| M   | NoColitis |
        | NC3    | 70.5| F   | NoColitis |
        | NC4    | 63.7| M   | NoColitis |
        | NC5    | 58.5| M   | NoColitis |
        | NC6    | 49.3| F   | NoColitis |
        | CT1    | 21.4| F   | Control   |
        | CT2    | 61.7| M   | Control   |
        | CT3    | 50.5| M   | Control   |
        | CT4    | 43.4| M   | Control   |
        | CT5    | 70.6| F   | Control   |
        | CT6    | 44.3| M   | Control   |
        | CT7    | 50.2| M   | Control   |
        | CT8    | 61.5| F   | Control   |

        ### Count the number of samples per Diagnosis

        ```toml
        [SampleInfo.envs.stats."N_Samples_per_Diagnosis (pie)"]
        plot_type = "pie"
        x = "sample"
        split_by = "Diagnosis"
        ```

        ![Samples_Diagnosis]({{output_baseurl}}/sampleinfo/SampleInfo/N_Samples_per_Diagnosis-pie-.png)

        What if we want a bar plot instead of a pie chart?

        ```toml
        [SampleInfo.envs.stats."N_Samples_per_Diagnosis (bar)"]
        plot_type = "bar"
        x = "Sample"
        split_by = "Diagnosis"
        ```

        ![Samples_Diagnosis_bar]({{output_baseurl}}/sampleinfo/SampleInfo/N_Samples_per_Diagnosis-bar-.png)

        ### Explore Age distribution

        The distribution of Age of all samples

        ```toml
        [SampleInfo.envs.stats."Age_distribution (histogram)"]
        plot_type = "histogram"
        x = "Age"
        ```

        ![Age_distribution]({{output_baseurl}}/sampleinfo/SampleInfo/Age_distribution-Histogram-.png)

        How about the distribution of Age in each Diagnosis, and make it
        violin + boxplot?

        ```toml
        [SampleInfo.envs.stats."Age_distribution_per_Diagnosis (violin + boxplot)"]
        y = "Age"
        x = "Diagnosis"
        plot_type = "violin"
        add_box = true
        ```

        ![Age_distribution_per_Diagnosis]({{output_baseurl}}/sampleinfo/SampleInfo/Age_distribution_per_Diagnosis-violin-boxplot-.png)

        How about Age distribution per Sex in each Diagnosis?

        ```toml
        [SampleInfo.envs.stats."Age_distribution_per_Sex_in_each_Diagnosis (boxplot)"]
        y = "Age"
        x = "Sex"
        split_by = "Diagnosis"
        plot_type = "box"
        ncol = 3
        devpars = {height = 450}
        ```

        ![Age_distribution_per_Sex_in_each_Diagnosis]({{output_baseurl}}/sampleinfo/SampleInfo/Age_distribution_per_Sex_in_each_Diagnosis-boxplot-.png)

    Input:
        infile%(required)s: {{Input.infile.help | indent: 8}}.
            **Required when [`LoadingRNAFromSeurat`](LoadingRNAFromSeurat.md) is not
            used in the pipeline.**
            **This is optional if [`LoadingRNAFromSeurat`](LoadingRNAFromSeurat.md) is
            used in the pipeline and no VDJ data is provided.**
            The input file should have the following columns.
            * Sample: A unique id for each sample.
            * TCRData/BCRData: The directory for single-cell TCR/BCR data for this
                sample.
                Specifically, it should contain filtered_contig_annotations.csv
                or all_contig_annotations.csv from cellranger.
            * RNAData: The directory for single-cell RNA data for this sample.
                Specifically, it should be able to be read by
                [`Seurat::Read10X()`](https://satijalab.org/seurat/reference/read10x) or
                [`Seurat::Read10X_h5()`](https://satijalab.org/seurat/reference/read10x_h5) or
                [`SeuratDisk::LoadLoom()`](https://rdrr.io/github/mojaveazure/seurat-disk/man/LoadLoom.html).
                See also <https://satijalab.org/seurat/reference/read10x>.
            * Other columns are optional and will be treated as metadata for
                each sample.
    """ % {  # noqa: E501
        "required": " (required)" if "LoadingRNAFromSeurat" not in config else ""
    }

    envs = {"exclude_cols": "TCRData,BCRData,RNAData"}


# When SampleInfo is used, it should always have VDJ data to be loaded
# if has_vdj is True
@when(SampleInfo and config.has_vdj, requires=SampleInfo)
@annotate.format_doc()
class ScRepLoading(ScRepLoading_):
    pass


VDJInput = ScRepLoading


# Input from Seurat object, it doesn't require SampleInfo
@when("LoadingRNAFromSeurat" in config)
class LoadingRNAFromSeurat(Proc):
    """Load RNA data from a Seurat object, instead of RNAData from SampleInfo

    Input:
        infile: An [RDS](https://rdrr.io/r/base/readRDS.html) or [qs/qs2](https://github.com/qsbase/qs2)
        format file containing a Seurat object.

    Envs:
        prepared (flag): Whether the Seurat object is well-prepared for the
            pipeline (so that SeuratPreparing process is not needed).
        clustered (flag): Whether the Seurat object is clustered, so that
            `SeuratClustering` (`SeuratClusteringOfAllCells`) process or
            `SeuratMap2Ref` is not needed.
            Force `prepared` to be `True` if this is `True`.
        sample: The column name in the metadata of the Seurat object that
            indicates the sample name.

    SeeAlso:
        - [Preparing the input](../preparing-input.md#single-cell-rna-seq-scrna-seq-data).
        - [Routes of the pipeline](../introduction.md#routes-of-the-pipeline).
    """  # noqa: E501

    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = biopipen_config.lang.rscript
    envs = {
        "prepared": False,
        "clustered": False,
        "sample": "Sample",
    }
    script = "file://scripts/LoadingRNAFromSeurat.R"


# Ensured by validate_config that either loads
RNAInput = LoadingRNAFromSeurat or SampleInfo


@when(
    (
        # Even when we load RNA-seq data from Seurat, we may still need SeuratPreparing
        # for QC, transformation, etc.
        "LoadingRNAFromSeurat" in config
        and not config.LoadingRNAFromSeurat.envs.prepared
    )
    or (
        # Or when we load RNA-seq data from SampleInfo
        SampleInfo
        and not LoadingRNAFromSeurat
    ),
    requires=RNAInput,
)
@annotate.format_doc()
class SeuratPreparing(SeuratPreparing_):
    """{{Summary}}

    See also [Preparing the input](../preparing-input.md#single-cell-rna-seq-scrna-seq-data).

    Metadata:
        Here is the demonstration of basic metadata for the `Seurat` object. Future
        processes will use it and/or add more metadata to the `Seurat` object.

        ![SeuratPreparing-metadata](images/SeuratPreparing-metadata.png)
    """  # noqa: E501


RNAInput = SeuratPreparing or RNAInput


# No matter "SeuratClusteringOfAllCells" is in the config or not
# if TOrBCellSelection is used, meaning input RNA data has T/B cells and non-T/B cells
@when("TOrBCellSelection" in config, requires=RNAInput)
@annotate.format_doc()
class SeuratClusteringOfAllCells(SeuratClustering_):
    """Cluster all cells, including T cells/non-T cells and B cells/non-Bcells
    using Seurat.

    This process will perform clustering on all cells using
    [`Seurat`](https://satijalab.org/seurat/) package.
    The clusters will then be used to select T/B cells by
    [`TOrBCellSelection`](TOrBCellSelection.md) process.

    {{*Summary.long}}

    /// Note
    If all your cells are all T/B cells ([`TOrBCellSelection`](TOrBCellSelection.md)
    is not set in configuration), you should not use this process.
    Instead, you should use [`SeuratClustering`](./SeuratClustering.md) process
    for unsupervised clustering, or [`SeuratMap2Ref`](./SeuratMap2Ref.md) process
    for supervised clustering.
    ///

    SeeAlso:
        - [SeuratClustering](./SeuratClustering.md)
    """


RNAInput = SeuratClusteringOfAllCells or RNAInput


@when(SeuratClusteringOfAllCells, requires=RNAInput)
@annotate.format_doc()
class ClusterMarkersOfAllCells(MarkersFinder_):
    """Markers for clusters of all cells.

    SeeAlso:
        - [ClusterMarkers](./ClusterMarkers.md)
        - [MarkersFinder](./MarkersFinder.md)
        - [biopipen.ns.scrna.MarkersFinder](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder)

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        ident_1 (hidden;readonly): {{Envs["ident_1"].help | indent: 12}}.
        ident_2 (hidden;readonly): {{Envs["ident_2"].help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
    """  # noqa: E501

    envs = {
        "cases": {"Cluster": {"group_by": "seurat_clusters"}},
        "marker_plots_defaults": {"order_by": "desc(avg_log2FC)"},
        "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
        "allmarker_plots": {"Top 10 markers of all clusters": {"plot_type": "heatmap"}},
    }
    order = 2


@when(
    SeuratClusteringOfAllCells and "TopExpressingGenesOfAllCells" in config,
    requires=RNAInput,
)
@annotate.format_doc()
class TopExpressingGenesOfAllCells(TopExpressingGenes_):
    """Top expressing genes for clusters of all cells.

    {{*Summary.long}}

    SeeAlso:
        - [TopExpressingGenes](./TopExpressingGenes.md)
        - [ClusterMarkers](./ClusterMarkers.md) for examples of enrichment plots

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        group_by (hidden;readonly): {{Envs["group_by"].help | indent: 12}}.
        ident (hidden;readonly): {{Envs.ident.help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
    """

    envs = {"cases": {"Cluster": {}}}
    order = 3


@when(
    "TOrBCellSelection" in config,
    requires=[RNAInput, VDJInput] if VDJInput else RNAInput,
)
@annotate.format_doc()
class TOrBCellSelection(TOrBCellSelection_):
    pass


RNAInput = TOrBCellSelection or RNAInput


@when("ModuleScoreCalculator" in config, requires=RNAInput)
@annotate.format_doc()
class ModuleScoreCalculator(ModuleScoreCalculator_):
    """{{Summary}}

    Metadata:
        The metadata of the `Seurat` object will be updated with the module scores:

        ![ModuleScoreCalculator-metadata](images/ModuleScoreCalculator-metadata.png)
    """  # noqa: E501

    input_data = lambda ch1: ch1.iloc[:, [0]]


RNAInput = ModuleScoreCalculator or RNAInput


@when("SeuratMap2Ref" in config, requires=RNAInput)
@annotate.format_doc()
class SeuratMap2Ref(SeuratMap2Ref_):
    """{{Summary}}

    Metadata:
        The metadata of the `Seurat` object will be updated with the cluster
        assignments (column name determined by `envs.name`):

        ![SeuratMap2Ref-metadata](images/SeuratClustering-metadata.png)
    """

    input_data = lambda ch1: ch1.iloc[:, [0]]


RNAInput = SeuratMap2Ref or RNAInput


@when(
    "SeuratClustering" in config
    or (
        "SeuratMap2Ref" not in config
        and "CellTypeAnnotation" not in config
        and (
            "LoadingRNAFromSeurat" not in config
            or not config.LoadingRNAFromSeurat.envs.clustered
        )
    ),
    requires=RNAInput,
)
@annotate.format_doc()
class SeuratClustering(SeuratClustering_):
    """Cluster all cells or selected T/B cells selected by `TOrBCellSelection`.

    If `[TOrBCellSelection]` is not set in the configuration, meaning
    all cells are T/B cells, this process will be run on all T/B cells. Otherwise,
    this process will be run on the selected T/B cells by
    [`TOrBCellSelection`](./TOrBCellSelection.md).

    SeeAlso:
        - [SeuratClusteringOfAllCells](./SeuratClusteringOfAllCells.md)

    Metadata:
        The metadata of the `Seurat` object will be updated with the cluster
        assignments:

        ![SeuratClustering-metadata](images/SeuratClustering-metadata.png)
    """

    input_data = lambda ch1: ch1.iloc[:, [0]]


RNAInput = SeuratClustering or RNAInput


@when("CellTypeAnnotation" in config, requires=RNAInput)
@annotate.format_doc()
class CellTypeAnnotation(CellTypeAnnotation_):
    """Annotate all or selected T/B cell clusters.

    {{*Summary}}

    The `<workdir>` is typically `./.pipen` and the `<pipline_name>` is `Immunopipe`
    by default.

    /// Note
    When supervised clustering [`SeuratMap2Ref`](./SeuratMap2Ref.md) is used, this
    process will be ignored.
    ///

    /// Note
    When cell types are annotated, the old `seurat_clusters` column will be renamed
    to `seurat_clusters_id`, and the new `seurat_clusters` column will be added.
    ///

    Metadata:
        When `envs.tool` is `direct` and `envs.cell_types` is empty, the metadata of
        the `Seurat` object will be kept as is.

        When `envs.newcol` is specified, the original `seurat_clusters` column will
        be kept is, and the annotated cell types will be saved in the new column.
        Otherwise, the original `seurat_clusters` column will be replaced by the
        annotated cell types and the original `seurat_clusters` column will be
        saved at `seurat_clusters_id`.

        ![CellTypeAnnotation-metadata](images/CellTypeAnnotation-metadata.png)

    """  # noqa: E501

    # Change the default to direct, which doesn't do any annotation
    envs = {"tool": "direct", "sctype_db": None}


RNAInput = CellTypeAnnotation or RNAInput


@when("SeuratSubClustering" in config, requires=RNAInput)
@annotate.format_doc()
class SeuratSubClustering(SeuratSubClustering_):
    """Sub-clustering for all or selected T/B cells.

    {{*Summary}}

    Metadata:
        The metadata of the `Seurat` object will be updated with the sub-clusters
        specified by names (keys) of `envs.cases`:

        ![SeuratSubClustering-metadata](images/SeuratSubClustering-metadata.png)
    """

    input_data = lambda ch1: ch1.iloc[:, [0]]


RNAInput = SeuratSubClustering or RNAInput


@annotate.format_doc(vars={"output_baseurl": TEST_OUTPUT_BASEURL})
class ClusterMarkers(MarkersFinder_):
    """Markers for clusters of all or selected T/B cells.

    This process is extended from [`MarkersFinder`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder)
    from the [`biopipen`](https://pwwang.github.io/biopipen) package.
    `MarkersFinder` is a `pipen` process that wraps the
    [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers)
    function, and performs enrichment analysis for the markers found.

    The enrichment analysis is done by [`enrichr`](https://maayanlab.cloud/Enrichr/).

    /// Note
    Since this process is extended from `MarkersFinder`, other environment variables from `MarkersFinder` are also available.
    However, they should not be used in this process. Other environment variables are used for more complicated cases for marker finding
    (See [`MarkersFinder`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder) for more details).

    If you are using `pipen-board` to run the pipeline
    (see [here](../running.md#run-the-pipeline-via-pipen-board) and
    [here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)),
    you may see the other environment variables of this process are hidden and readonly.
    ///

    SeeAlso:
        - [MarkersFinder](./MarkersFinder.md)
        - [ClusterMarkersOfAllCells](./ClusterMarkersOfAllCells.md)
        - [biopipen.ns.scrna.MarkersFinder](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder)

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        group_by (hidden;readonly): {{Envs["group_by"].help | indent: 12}}.
        ident_1 (hidden;readonly): {{Envs["ident_1"].help | indent: 12}}.
        ident_2 (hidden;readonly): {{Envs["ident_2"].help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.

    Examples:
        ### Visualize Log2 Fold Change of Markers

        ```toml
        [ClusterMarkers.envs.marker_plots."Volcano Plot (log2FC)"]
        plot_type = "volcano_log2fc"
        ```

        ![Volcano Plot (log2FC)]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Volcano-Plot-log2FC-.png)

        ### Visualize differential percentage of expression of Markers

        ```toml
        [ClusterMarkers.envs.marker_plots."Volcano Plot (pct_diff)"]
        plot_type = "volcano_pct"
        ```

        ![Volcano Plot (pct_diff)]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Volcano-Plot-diff_pct-.png)

        ### Visualize Average Expression of Markers with Dot Plot

        ```toml
        [ClusterMarkers.envs.marker_plots."Dot Plot (AvgExp)"]
        plot_type = "dotplot"
        order_by = "desc(avg_log2FC)"
        ```

        ![Dot Plot (AvgExp)]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Dot-Plot.png)

        ### Visualize Average Expression of Markers with Heatmap

        ```toml
        [ClusterMarkers.envs.marker_plots."Heatmap (AvgExp)"]
        plot_type = "heatmap"
        order_by = "desc(avg_log2FC)"
        ```

        ![Heatmap (AvgExp)]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Heatmap-of-Expressions-of-Top-Markers.png)

        ### Visualize Expression of Markers with Violin Plots

        ```toml
        [ClusterMarkers.envs.marker_plots."Violin Plots"]
        plot_type = "violin"
        ```

        ![Violin Plots]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Violin-Plots-for-Top-Markers.png)

        ### Visualize enrichment analysis results with Bar/EnrichMap/Network/WordCloud Plots

        ```toml
        # Visualize enrichment of markers
        [ClusterMarkers.envs.enrich_plots."Bar Plot"]  # Default
        plot_type = "bar"

        [ClusterMarkers.envs.enrich_plots."Network"]
        plot_type = "network"

        [ClusterMarkers.envs.enrich_plots."Enrichmap"]
        plot_type = "enrichmap"

        [ClusterMarkers.envs.enrich_plots."Word Cloud"]
        plot_type = "wordcloud"
        ```

        ![Bar Plot]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Bar-Plot.png)
        ![Network]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Network.png)
        ![Enrichmap]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Enrichmap.png)
        ![Word Cloud]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/enrich.MSigDB_Hallmark_2020.Word-Cloud.png)

        ### Visualize top markers of all clusters with Heatmap

        ```toml
        [ClusterMarkers.envs.allmarker_plots."Top 10 markers of all clusters"]
        plot_type = "heatmap"
        ```

        ![Top 10 markers of all clusters]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/Top-10-markers-of-all-clusters.png)

        ### Visualize Log2 Fold Change of all markers

        ```toml
        [ClusterMarkers.envs.allmarker_plots."Log2 Fold Change of all markers"]
        plot_type = "heatmap_log2fc"
        subset_by = "seurat_clusters"
        ```

        ![Log2 Fold Change of all markers]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/Log2FC-of-all-clusters.png)

        ### Visualize all markers in all clusters with Jitter Plots

        ```toml
        [ClusterMarkers.envs.allmarker_plots."Jitter Plots of all markers"]
        plot_type = "jitter"
        subset_by = "seurat_clusters"
        ```

        ![Jitter Plots of all markers]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/Jitter-Plots-for-all-clusters.png)

        ### Visualize all enrichment analysis results of all clusters

        ```toml
        [ClusterMarkers.envs.allenrich_plots."Heatmap of enriched terms of all clusters"]
        plot_type = "heatmap"
        ```

        ![Heatmap of enriched terms of all clusters]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Enrichments-/allenrich.MSigDB_Hallmark_2020.Heatmap-of-enriched-terms-of-all-clusters.png)

        ### Overlapping markers

        ```toml
        [ClusterMarkers.envs.overlaps."Overlapping Markers"]
        plot_type = "venn"
        ```

        ![Overlapping Markers]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-Overlaps-/Overlapping-Markers.png)

    """  # noqa: E501

    requires = RNAInput
    envs = {
        "cases": {"Cluster": {"group_by": "seurat_clusters"}},
        "marker_plots_defaults": {"order_by": "desc(avg_log2FC)"},
        "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
        "allmarker_plots": {"Top 10 markers of all clusters": {"plot_type": "heatmap"}},
    }
    order = 2


@when("TopExpressingGenes" in config, requires=RNAInput)
@annotate.format_doc()
class TopExpressingGenes(TopExpressingGenes_):
    """Top expressing genes for clusters of all or selected T/B cells.

    {{*Summary.long}}

    This process finds the top expressing genes of clusters of T/B cells, and also
    performs the enrichment analysis against the genes.

    The enrichment analysis is done by
    [`enrichr`](https://maayanlab.cloud/Enrichr/).

    /// Note
    There are other environment variables also available. However, they should not
    be used in this process. Other environment variables are used for more
    complicated cases for investigating top genes
    (See [`biopipen.ns.scrna.TopExpressingGenes`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.TopExpressingGenes) for more details).

    If you are using `pipen-board` to run the pipeline
    (see [here](../running.md#run-the-pipeline-via-pipen-board) and
    [here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)),
    you may see the other environment variables of this process are hidden and
    readonly.
    ///

    SeeAlso:
        - [TopExpressingGenesOfAllCells](./TopExpressingGenesOfAllCells.md)
        - [ClusterMarkers](./ClusterMarkers.md) for examples of enrichment plots

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        group_by (hidden;readonly): {{Envs["group_by"].help | indent: 12}}.
        ident (hidden;readonly): {{Envs.ident.help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
    """  # noqa: E501

    envs = {"cases": {"Cluster": {}}}
    order = 3


@when(VDJInput, requires=[VDJInput, RNAInput])
class ScRepCombiningExpression(ScRepCombiningExpression_):
    pass


CombinedInput = ScRepCombiningExpression or RNAInput


@when(VDJInput and "TCRClustering" in config, requires=CombinedInput)
@annotate.format_doc()
class TCRClustering(TCRClustering_):
    """{{Summary.short}}

    You can disable this by remving the whole sections of
    TCRClustering in the config file.

    {{*Summary.long}}
    """

    input_data = lambda ch1: ch1.iloc[:, [0]]
    order = 4


CombinedInput = TCRClustering or CombinedInput


@when(VDJInput and "TESSA" in config, requires=CombinedInput)
@annotate.format_doc()
class TESSA(TESSA_):
    """{{Summary}}

    /// Note
    The dependencies of TESSA are not included in the docker image of immunopipe
    with tag without `-full` suffix. If you want to use TESSA, please use the
    docker image with tag with `-full` suffix, or install the dependencies manually.
    ///

    Metadata:
        The metadata of the `Seurat` object will be updated with the TESSA clusters
        and the cluster sizes:

        ![TESSA-metadata](images/TESSA-metadata.png)
    """

    order = 5


CombinedInput = TESSA or CombinedInput


@when(
    "CellCellCommunication" in config or "CellCellCommunicationPlots" in config,
    requires=CombinedInput,
)
class CellCellCommunication(CellCellCommunication_):
    order = 7


@when(
    "CellCellCommunication" in config or "CellCellCommunicationPlots" in config,
    requires=CellCellCommunication,
)
class CellCellCommunicationPlots(CellCellCommunicationPlots_):
    pass


class SeuratClusterStats(SeuratClusterStats_):
    requires = CombinedInput
    order = -1
    envs_depth = 3
    envs = {
        "dimplots": {
            "Dimensional reduction plot": {
                "label": True,
            },
        },
    }
    if VDJInput:
        envs["dimplots"]["VDJ Presence"] = {
            "group_by": "VDJ_Presence",
        }


@when(VDJInput, requires=CombinedInput)
class ClonalStats(ClonalStats_):
    envs_depth = 3
    order = 8


@when("ScFGSEA" in config, requires=CombinedInput)
class ScFGSEA(ScFGSEA_):
    order = 9


@when("PseudoBulkDEG" in config, requires=CombinedInput)
@annotate.format_doc()
class PseudoBulkDEG(PseudoBulkDEG_):
    """{{Summary}}

    SeeAlso:
        - [biopipen.ns.scrna.PseudoBulkDEG](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.PseudoBulkDEG)
        - [ClusterMarkers](./ClusterMarkers.md) for examples of marker and enrichment plots
    """  # noqa: E501
    order = 10


@when("MarkersFinder" in config, requires=CombinedInput)
@annotate.format_doc()
class MarkersFinder(MarkersFinder_):
    """{{Summary.short}}

    `MarkersFinder` is a process that wraps the
    [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers)
    function, and performs enrichment analysis for the markers found.

    SeeAlso:
        - [biopipen.ns.scrna.MarkersFinder](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.MarkersFinder)
        - [ClusterMarkers](./ClusterMarkers.md)

    Envs:
        mutaters: {{Envs.mutaters.help | indent: 12}}.
            See also
            [mutating the metadata](../configurations.md#mutating-the-metadata).

    Examples:
        The examples are for more general use of `MarkersFinder`, in order to
        demonstrate how the final cases are constructed.

        Suppose we have a metadata like this:

        | id | seurat_clusters | Group |
        |----|-----------------|-------|
        | 1  | 1               | A     |
        | 2  | 1               | A     |
        | 3  | 2               | A     |
        | 4  | 2               | A     |
        | 5  | 3               | B     |
        | 6  | 3               | B     |
        | 7  | 4               | B     |
        | 8  | 4               | B     |

        ### Default

        By default, `group_by` is `seurat_clusters`, and `ident_1` and `ident_2`
        are not specified. So markers will be found for all clusters in the manner
        of "cluster vs rest" comparison.

        - Cluster
            - 1 (vs 2, 3, 4)
            - 2 (vs 1, 3, 4)
            - 3 (vs 1, 2, 4)
            - 4 (vs 1, 2, 3)

        Each case will have the markers and the enrichment analysis for the
        markers as the results.

        ### With `each` group

        `each` is used to separate the cells into different cases. `group_by`
        is still `seurat_clusters`.

        ```toml
        [<Proc>.envs]
        group_by = "seurat_clusters"
        each = "Group"
        ```

        - A:Cluster
            - 1 (vs 2)
            - 2 (vs 1)
        - B:Cluster
            - 3 (vs 4)
            - 4 (vs 3)

        ### With `ident_1` only

        `ident_1` is used to specify the first group of cells to compare.
        Then the rest of the cells in the case are used for `ident_2`.

        ```toml
        [<Proc>.envs]
        group_by = "seurat_clusters"
        ident_1 = "1"
        ```

        - Cluster
            - 1 (vs 2, 3, 4)

        ### With both `ident_1` and `ident_2`

        `ident_1` and `ident_2` are used to specify the two groups of cells to
        compare.

        ```toml
        [<Proc>.envs]
        group_by = "seurat_clusters"
        ident_1 = "1"
        ident_2 = "2"
        ```

        - Cluster
            - 1 (vs 2)

        ### Multiple cases

        ```toml
        [<Proc>.envs.cases]
        c1_vs_c2 = {ident_1 = "1", ident_2 = "2"}
        c3_vs_c4 = {ident_1 = "3", ident_2 = "4"}
        ```

        - DEFAULT:c1_vs_c2
            - 1 (vs 2)
        - DEFAULT:c3_vs_c4
            - 3 (vs 4)

        The `DEFAULT` section name will be ignored in the report. You can specify
        a section name other than `DEFAULT` for each case to group them
        in the report.
    """  # noqa: E501

    order = 11


@when(VDJInput and "CDR3AAPhyschem" in config, requires=CombinedInput)
class CDR3AAPhyschem(CDR3AAPhyschem_):
    order = 12


if "ScrnaMetabolicLandscape" in config or just_loading:
    anno = annotate(ScrnaMetabolicLandscape)
    anno.Args.metafile.attrs["readonly"] = True
    anno.Args.metafile.attrs["hidden"] = True
    anno.Args.is_seurat.attrs["readonly"] = True
    anno.Args.is_seurat.attrs["hidden"] = True
    anno.Args.is_seurat.attrs["flag"] = True
    anno.Args.is_seurat.attrs["default"] = True
    anno.Args.is_seurat.attrs["value"] = True

    if just_loading:
        scrna_metabolic_landscape = ScrnaMetabolicLandscape(
            is_seurat=True,
            noimpute=False,
        )
    else:
        scrna_metabolic_landscape = ScrnaMetabolicLandscape(is_seurat=True)

    scrna_metabolic_landscape.p_input.requires = CombinedInput
    scrna_metabolic_landscape.p_input.order = 99
