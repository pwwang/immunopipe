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
    CDR3Clustering as CDR3Clustering_,
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
    Slingshot as Slingshot_,
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

    def decorator(cls: Type[Proc]) -> Type[Proc] | None:
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

        ![Samples_Diagnosis]({{output_baseurl}}/sampleinfo/SampleInfo/N_Samples_per_Diagnosis-pie.png)

        What if we want a bar plot instead of a pie chart?

        ```toml
        [SampleInfo.envs.stats."N_Samples_per_Diagnosis (bar)"]
        plot_type = "bar"
        x = "Sample"
        split_by = "Diagnosis"
        ```

        ![Samples_Diagnosis_bar]({{output_baseurl}}/sampleinfo/SampleInfo/N_Samples_per_Diagnosis-bar.png)

        ### Explore Age distribution

        The distribution of Age of all samples

        ```toml
        [SampleInfo.envs.stats."Age_distribution (histogram)"]
        plot_type = "histogram"
        x = "Age"
        ```

        ![Age_distribution]({{output_baseurl}}/sampleinfo/SampleInfo/Age_distribution-Histogram.png)

        How about the distribution of Age in each Diagnosis, and make it
        violin + boxplot?

        ```toml
        [SampleInfo.envs.stats."Age_distribution_per_Diagnosis (violin + boxplot)"]
        y = "Age"
        x = "Diagnosis"
        plot_type = "violin"
        add_box = true
        ```

        ![Age_distribution_per_Diagnosis]({{output_baseurl}}/sampleinfo/SampleInfo/Age_distribution_per_Diagnosis-violin-boxplot.png)

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

        ![Age_distribution_per_Sex_in_each_Diagnosis]({{output_baseurl}}/sampleinfo/SampleInfo/Age_distribution_per_Sex_in_each_Diagnosis-boxplot.png)

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

    Description:
        This is the entrance of the pipeline: it lists the sample information
        given in the input file and performs descriptive statistics and plots on
        it (`envs.stats`). The statistics and plots are produced with `dplyr` and
        `plotthis` (R); no upstream analysis tool is wrapped.

    Base class:
        `biopipen.ns.delim.SampleInfo`

    Deviations:
        `exclude_cols` is overridden: it is `None` in the base and is set to
        `TCRData,BCRData,RNAData` here, so the three data-path columns are kept
        out of the statistics and out of the report.
    """ % {  # noqa: E501
        "required": " (required)" if "LoadingRNAFromSeurat" not in config else ""
    }

    envs = {"exclude_cols": "TCRData,BCRData,RNAData"}


# When SampleInfo is used, it should always have VDJ data to be loaded
# if has_vdj is True
@when(SampleInfo and config.has_vdj, requires=SampleInfo)  # type: ignore
@annotate.format_doc()
class ScRepLoading(ScRepLoading_):
    """Load the single cell TCR/BCR data into a `scRepertoire` compatible object

    Description:
        Loads the TCR/BCR (VDJ) data of each sample into a `scRepertoire`
        compatible object, so that `ScRepCombiningExpression` can later combine
        the repertoire data with the expression data. The loading is done by
        `scRepertoire` (R), through `scRepertoire::loadContigs()`.

    Base class:
        `biopipen.ns.tcr.ScRepLoading`

    Deviations:
        The parameters are passed through to `scRepertoire` unchanged: this
        process overrides no inherited parameter and adds none.
    """


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
            Multiple columns will be concatenated with `_` to form the sample name.
        mutaters (type=json): The mutaters to mutate the metadata
            Keys are the names of the mutaters and values are the R expressions
            passed by `dplyr::mutate()` to mutate the metadata.
        subset: An expression to subset the cells, will be passed to `dplyr::filter()`.
            This will be applied after mutating the metadata.
        ncores (type=int): The number of threads used to load/save the Seurat object.

    SeeAlso:
        - [Preparing the input](../preparing-input.md#single-cell-rna-seq-scrna-seq-data).
        - [Routes of the pipeline](../introduction.md#routes-of-the-pipeline).

    Description:
        Loads the RNA data from a pre-existing `Seurat` object (an RDS or
        qs/qs2 file) instead of the `RNAData` directories listed by
        `SampleInfo`. This is not a wrapper of an upstream analysis tool: the
        process is immunopipe's own and runs an R script that ships with
        immunopipe (`immunopipe/scripts/LoadingRNAFromSeurat.R`), which reads
        and writes the object with `tidyseurat` and `qs2`.

    Base class:
        `biopipen.core.proc.Proc` - biopipen's bare process class, which
        declares no parameters of its own.

    Deviations:
        All six options of this process are added by immunopipe, since the base
        class declares none: `prepared` (default `False`), `clustered`
        (default `False`), `sample` (default `Sample`), `mutaters`
        (default `{}`), `subset` (default `None`) and `ncores` (default `1`).
    """  # noqa: E501

    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = biopipen_config.lang.rscript
    envs = {
        "prepared": False,
        "clustered": False,
        "sample": "Sample",
        "mutaters": {},
        "subset": None,
        "ncores": biopipen_config.misc.ncores,
    }
    script = "file://scripts/LoadingRNAFromSeurat.R"


# Ensured by validate_config that either loads
RNAInput = LoadingRNAFromSeurat or SampleInfo


@when(
    (
        # Even when we load RNA-seq data from Seurat, we may still need SeuratPreparing
        # for QC, transformation, etc.
        "LoadingRNAFromSeurat" in config
        and not config.LoadingRNAFromSeurat.envs.prepared  # type: ignore
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

    Description:
        Loads the scRNA-seq data, prepares it (normalization and integration of
        the samples) and applies quality control to it, using `Seurat` (R) -
        `CreateSeuratObject()` and `Read10X()` for the loading, and the cell-
        and gene-level filters given under `envs.cell_qc` and `envs.gene_qc`.

    Base class:
        `biopipen.ns.scrna.SeuratPreparing`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
    """  # noqa: E501

    # Don't export the RDS/qs file
    export = False


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

    Description:
        Clusters all the cells of the object - T cells together with non-T
        cells, or B cells together with non-B cells - so that
        `TOrBCellSelection` has a clustering of the whole dataset to select the
        T/B cells from. The clustering is done by `Seurat` (R), with
        `FindNeighbors()`, `FindClusters()` and `RunUMAP()`.

    Base class:
        `biopipen.ns.scrna.SeuratClustering`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
    """


RNAInput = SeuratClusteringOfAllCells or RNAInput


@when(SeuratClusteringOfAllCells, requires=RNAInput)  # type: ignore
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

    Description:
        Finds the marker genes of every cluster found by
        `SeuratClusteringOfAllCells` and runs an enrichment analysis on them.
        Markers are found by `Seurat::FindMarkers()` and the enrichment is done
        by `enrichr`.

    Base class:
        `biopipen.ns.scrna.MarkersFinder`

    Deviations:
        Four inherited parameters are changed. `sigmarkers` is set to
        `p_val_adj < 0.05 & avg_log2FC > 0` instead of `p_val_adj < 0.05`, so
        that only up-regulated markers are kept. `allmarker_plots` is set to a
        `heatmap` of the top 10 markers of all clusters.
        `marker_plots_defaults` gains `order_by = "desc(avg_log2FC)"`. `cases`
        is empty in the base, so the base runs no marker-finding case by
        default; immunopipe defines the `Cluster` case with `group_by = None`.
    """  # noqa: E501

    envs = {
        "cases": {"Cluster": {"group_by": None}},
        "marker_plots_defaults": {"order_by": "desc(avg_log2FC)"},
        "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
        "allmarker_plots": {"Top 10 markers of all clusters": {"plot_type": "heatmap"}},
    }
    order = 2


@when(
    SeuratClusteringOfAllCells  # type: ignore
    and "TopExpressingGenesOfAllCells" in config,
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

    Description:
        Finds the top expressing genes of every cluster found by
        `SeuratClusteringOfAllCells` and runs an enrichment analysis on them.
        The top expressing genes are computed by `Seurat` and the enrichment is
        done by `enrichr`.

    Base class:
        `biopipen.ns.scrna.TopExpressingGenes`

    Deviations:
        `cases` is empty in the base, so the base computes no case by default;
        immunopipe defines the `Cluster` case explicitly.
    """

    envs = {"cases": {"Cluster": {}}}
    order = 3


@when(
    "TOrBCellSelection" in config,
    requires=[RNAInput, VDJInput] if VDJInput else RNAInput,  # type: ignore
)
@annotate.format_doc()
class TOrBCellSelection(TOrBCellSelection_):
    """Separate T and non-T cells and select T cells; or separate B and
    non-B cells and select B cells.

    Description:
        Separates T from non-T cells (or B from non-B cells) and keeps the T/B
        cells for the downstream analysis. The selection uses the expression
        values of `envs.indicator_genes` and, unless `envs.ignore_vdj` is set,
        the clonotype percentage of the clusters; when no `envs.selector` is
        given, `stats::kmeans` (R, with K=2) separates the two groups. This is
        not a wrapper of an upstream analysis tool: the process is immunopipe's
        own, extending an immunopipe class and running immunopipe's own R script
        (`immunopipe/scripts/TOrBCellSelection.R`).

    Base class:
        `immunopipe.inhouse.TOrBCellSelection` - an immunopipe class rather than
        a biopipen one.

    Deviations:
        The four parameters are immunopipe's own and are used as the base class
        declares them; this process overrides none of them and adds none.
    """


RNAInput = TOrBCellSelection or RNAInput


@when(
    "SeuratClustering" in config
    or "CellTypeAnnotation" in config
    or (
        "SeuratMap2Ref" not in config
        and (
            "LoadingRNAFromSeurat" not in config
            or not config.LoadingRNAFromSeurat.envs.clustered  # type: ignore
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

    /// Note

    If you have other annotation processes, including
    [`SeuratMap2Ref`](./SeuratMap2Ref.md) process or
    [`CellTypeAnnotation`](./CellTypeAnnotation.md) process enabled in the same run,
    you can specify a different name for the column to store the cluster information
    using `envs.ident`, so that the results from different
    annotation processes won't overwrite each other.

    ///

    SeeAlso:
        - [SeuratClusteringOfAllCells](./SeuratClusteringOfAllCells.md)

    Metadata:
        The metadata of the `Seurat` object will be updated with the cluster
        assignments:

        ![SeuratClustering-metadata](images/SeuratClustering-metadata.png)

    Description:
        Clusters the cells to be analysed - all cells, or the T/B cells selected
        by `TOrBCellSelection` - using `Seurat` (R), with `FindNeighbors()`,
        `FindClusters()` and `RunUMAP()`.

    Base class:
        `biopipen.ns.scrna.SeuratClustering`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
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

    If you have other annotation processes, including [`SeuratClustering`](./SeuratClustering.md)
    process or [`SeuratMap2Ref`](./SeuratMap2Ref.md) process enabled in the same run,
    you may want to specify a different name for the column to store the annotated cell types
    using `envs.anno_col`, so that the results from different annotation processes won't overwrite each other.

    ///

    /// Attention

    If you are running the pipeline with the Docker image, following tools are not available in the Docker image:

    - Direct assignment: `direct`, `cell`
    - Marker-based: `ScType`, `hitype`, `scSorter`, `SCINA`, `CelliD`, `UCell`, `AUCell`, `GSVA`, `singscore`, `SCSA`, `MACA`
    - Model-based: `celltypist`, `SingleR`
    - LLM-based: `LLMCelltype`, `mLLMCelltype`, `LICT`
    - Reference-based: `scmap`, `CHETAH`, `scClassify`, `MapQuery`

    ///

    Metadata:
        When `envs.tool` is `direct` and `envs.cell_types` is empty, the metadata of
        the `Seurat` object will be kept as is.

        ![CellTypeAnnotation-metadata](images/CellTypeAnnotation-metadata.png)

    Description:
        Annotates the cells or the clusters with the annotation backend chosen
        with `envs.tool`. `Seurat` holds the object; the annotation itself is
        done by the selected backend, such as `celltypist`, `SingleR` or
        `hitype`.

    Base class:
        `biopipen.ns.scrna.CellTypeAnnotation`

    Deviations:
        `tool` is overridden: the base defaults to `hitype`, immunopipe sets it
        to `direct`, which assigns cell types without running any annotation
        tool. `sctype_db` (default `None`) is added; the base only offers the
        nested `envs.sctype.db`. immunopipe also overrides `input_data`
        (`lambda ch1: ch1.iloc[:, [0]]`) so that the `Seurat` object is taken
        from the first input channel.
    """  # noqa: E501

    # Change the default to direct, which doesn't do any annotation
    envs = {"tool": "direct", "sctype_db": None}
    input_data = lambda ch1: ch1.iloc[:, [0]]


RNAInput = CellTypeAnnotation or RNAInput


@when("SeuratMap2Ref" in config, requires=RNAInput)
@annotate.format_doc()
class SeuratMap2Ref(SeuratMap2Ref_):
    """{{Summary}}

    /// Note

    If you have other annotation processes, including [`SeuratClustering`](./SeuratClustering.md)
    process or [`CellTypeAnnotation`](./CellTypeAnnotation.md) process enabled in the same run,
    you may want to specify a different name for the column to store the mapped cluster information
    using `envs.ident`, so that the results from different annotation processes won't overwrite each other.

    ///

    Metadata:
        The metadata of the `Seurat` object will be updated with the cluster
        assignments (column name determined by `envs.name`):

        ![SeuratMap2Ref-metadata](images/SeuratClustering-metadata.png)

    Description:
        Maps the object onto a reference `Seurat` object and transfers the
        reference labels to the query cells (supervised analysis), using
        `Seurat` (R): `FindTransferAnchors()` and `MapQuery()`.

    Base class:
        `biopipen.ns.scrna.SeuratMap2Ref`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
    """  # noqa: E501

    input_data = lambda ch1: ch1.iloc[:, [0]]
    # Don't export the RDS/qs file
    export = False


RNAInput = SeuratMap2Ref or RNAInput


@when("SeuratSubClustering" in config, requires=RNAInput)
@annotate.format_doc()
class SeuratSubClustering(SeuratSubClustering_):
    """Sub-clustering for all or selected T/B cells.

    {{*Summary}}

    Metadata:
        The metadata of the `Seurat` object will be updated with the sub-clusters
        specified by names (keys) of `envs.cases`:

        ![SeuratSubClustering-metadata](images/SeuratSubClustering-metadata.png)

    Description:
        Sub-clusters the selected cells or clusters, using `Seurat` (R),
        through `Seurat::FindSubCluster()`.

    Base class:
        `biopipen.ns.scrna.SeuratSubClustering`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
    """

    input_data = lambda ch1: ch1.iloc[:, [0]]


RNAInput = SeuratSubClustering or RNAInput


@when("Slingshot" in config, requires=RNAInput)
@annotate.format_doc()
class Slingshot(Slingshot_):
    """Trajectory inference using Slingshot

    Description:
        Infers the cell lineages and the pseudotime from the clustering, using
        the `slingshot` package (R/Bioconductor).

    Base class:
        `biopipen.ns.scrna.Slingshot`

    Deviations:
        `outtype` (default `qs2`) is added, so that the trajectories are written
        to a qs2 file; the base class has no such parameter.
    """

    envs = {"outtype": "qs2"}


RNAInput = Slingshot or RNAInput


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

        ![Volcano Plot (log2FC)]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Volcano-Plot-log2FC.png)

        ### Visualize differential percentage of expression of Markers

        ```toml
        [ClusterMarkers.envs.marker_plots."Volcano Plot (pct_diff)"]
        plot_type = "volcano_pct"
        ```

        ![Volcano Plot (pct_diff)]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-c1/markers.Volcano-Plot-diff_pct.png)

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

        ![Top 10 markers of all clusters]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers/Top-10-markers-of-all-clusters.png)

        ### Visualize Log2 Fold Change of all markers

        ```toml
        [ClusterMarkers.envs.allmarker_plots."Log2 Fold Change of all markers"]
        plot_type = "heatmap_log2fc"
        subset_by = "seurat_clusters"
        ```

        ![Log2 Fold Change of all markers]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers/Log2FC-of-all-clusters.png)

        ### Visualize all markers in all clusters with Jitter Plots

        ```toml
        [ClusterMarkers.envs.allmarker_plots."Jitter Plots of all markers"]
        plot_type = "jitter"
        subset_by = "seurat_clusters"
        ```

        ![Jitter Plots of all markers]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Markers/Jitter-Plots-for-all-clusters.png)

        ### Visualize all enrichment analysis results of all clusters

        ```toml
        [ClusterMarkers.envs.allenrich_plots."Heatmap of enriched terms of all clusters"]
        plot_type = "heatmap"
        ```

        ![Heatmap of enriched terms of all clusters]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-All-Enrichments/allenrich.MSigDB_Hallmark_2020.Heatmap-of-enriched-terms-of-all-clusters.png)

        ### Overlapping markers

        ```toml
        [ClusterMarkers.envs.overlaps."Overlapping Markers"]
        plot_type = "venn"
        ```

        ![Overlapping Markers]({{output_baseurl}}/clustermarkers/ClusterMarkers/sampleinfo.markers/Cluster/seurat_clusters-Overlaps/Overlapping-Markers.png)

    Description:
        Finds the marker genes of every cluster of the T/B cells (or of all
        cells) and runs an enrichment analysis on them. Markers are found by
        `Seurat::FindMarkers()` and the enrichment is done by `enrichr`.

    Base class:
        `biopipen.ns.scrna.MarkersFinder`

    Deviations:
        Four inherited parameters are changed. `sigmarkers` is set to
        `p_val_adj < 0.05 & avg_log2FC > 0` instead of `p_val_adj < 0.05`, so
        that only up-regulated markers are kept. `allmarker_plots` is set to a
        `heatmap_log2fc` of the top 5 markers of each cluster with
        `cutoff = 0.05`. `marker_plots_defaults` gains
        `order_by = "desc(avg_log2FC)"`. `cases` is empty in the base, so the
        base runs no marker-finding case by default; immunopipe defines the
        `Cluster` case with `group_by = None`.
    """  # noqa: E501

    requires = RNAInput  # type: ignore
    envs = {
        "cases": {"Cluster": {"group_by": None}},
        "marker_plots_defaults": {"order_by": "desc(avg_log2FC)"},
        "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
        "allmarker_plots": {
            "Top 5 markers of each cluster": {
                "plot_type": "heatmap_log2fc",
                "select": 5,
                "cutoff": 0.05,
            },
        },
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

    Description:
        Finds the top expressing genes of every cluster of the T/B cells (or of
        all cells) and runs an enrichment analysis on them. The top expressing
        genes are computed by `Seurat` and the enrichment is done by `enrichr`.

    Base class:
        `biopipen.ns.scrna.TopExpressingGenes`

    Deviations:
        `cases` is empty in the base, so the base computes no case by default;
        immunopipe defines the `Cluster` case explicitly.
    """  # noqa: E501

    envs = {"cases": {"Cluster": {}}}
    order = 3


@when("ModuleScoreCalculator" in config, requires=RNAInput)
@annotate.format_doc()
class ModuleScoreCalculator(ModuleScoreCalculator_):
    """{{Summary}}

    Metadata:
        The metadata of the `Seurat` object will be updated with the module scores:

        ![ModuleScoreCalculator-metadata](images/ModuleScoreCalculator-metadata.png)

    Description:
        Calculates the module scores of each cell, using `Seurat` (R), through
        `Seurat::AddModuleScore()` and `Seurat::CellCycleScoring()`.

    Base class:
        `biopipen.ns.scrna.ModuleScoreCalculator`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
    """  # noqa: E501


RNAInput = ModuleScoreCalculator or RNAInput


@when(VDJInput, requires=[VDJInput, RNAInput])  # type: ignore
class ScRepCombiningExpression(ScRepCombiningExpression_):
    """Combine the scTCR/BCR data with the expression data

    Description:
        Combines the repertoire data with the expression data, so that the
        clonotype information is available in the metadata of the `Seurat`
        object. It is done by `scRepertoire` (R), through
        `scRepertoire::combineExpression()`.

    Base class:
        `biopipen.ns.tcr.ScRepCombiningExpression`

    Deviations:
        The parameters are passed through to `scRepertoire` unchanged: this
        process overrides no inherited parameter and adds none.
    """


CombinedInput = ScRepCombiningExpression or RNAInput


@when(
    VDJInput and "CDR3Clustering" in config,  # type: ignore
    requires=CombinedInput,
)
@annotate.format_doc()
class CDR3Clustering(CDR3Clustering_):
    """Cluster the TCR/BCR clones by their CDR3 sequences

    Description:
        Clusters the TCR/BCR clones by the similarity of their CDR3 sequences,
        so that clones with similar receptors end up in the same cluster. The
        clustering is done by `ClusTCR` (Python) or by `GIANA`, whichever is
        selected with `envs.tool`.

    Base class:
        `biopipen.ns.tcr.CDR3Clustering`

    Deviations:
        The parameters are passed through to the selected tool unchanged: this
        process overrides no inherited parameter and adds none.
    """

    input_data = lambda ch1: ch1.iloc[:, [0]]
    order = 4


CombinedInput = CDR3Clustering or CombinedInput


@when(VDJInput and "TESSA" in config, requires=CombinedInput)  # type: ignore
@annotate.format_doc()
class TESSA(TESSA_):
    """{{Summary}}

    Metadata:
        The metadata of the `Seurat` object will be updated with the TESSA clusters
        and the cluster sizes:

        ![TESSA-metadata](images/TESSA-metadata.png)

    Description:
        Runs TESSA, a Bayesian model that integrates T cell receptor (TCR)
        sequence profiling with transcriptomes to find phenotype-associated TCR
        clusters. It is done by TESSA (Python), whose encoder and model are
        shipped with `biopipen`.

    Base class:
        `biopipen.ns.tcr.TESSA`

    Deviations:
        The parameters are passed through to TESSA unchanged: this process
        overrides no inherited parameter and adds none.
    """

    order = 5


CombinedInput = TESSA or CombinedInput


@when(
    "CellCellCommunication" in config or "CellCellCommunicationPlots" in config,
    requires=CombinedInput,
)
class CellCellCommunication(CellCellCommunication_):
    """Cell-cell communication inference

    Description:
        Infers cell-cell communication between the cell groups, based on the
        expression of ligand-receptor pairs. It is done by `LIANA` (Python),
        which offers a number of inference methods; `envs.method` selects the
        one to use and defaults to `cellchat`.

    Base class:
        `biopipen.ns.scrna.CellCellCommunication`

    Deviations:
        The parameters are passed through to `LIANA` unchanged: this process
        overrides no inherited parameter and adds none.
    """

    order = 7


@when(
    "CellCellCommunication" in config or "CellCellCommunicationPlots" in config,
    requires=CellCellCommunication,
)
class CellCellCommunicationPlots(CellCellCommunicationPlots_):
    """Visualization for cell-cell communication inference.

    Description:
        Draws the plots of the cell-cell communication results produced by
        `CellCellCommunication`, using `scplotter` (R), through
        `scplotter::CCCPlot()`.

    Base class:
        `biopipen.ns.scrna.CellCellCommunicationPlots`

    Deviations:
        The parameters are passed through to `scplotter` unchanged: this
        process overrides no inherited parameter and adds none.
    """


class SeuratClusterStats(SeuratClusterStats_):
    """Statistics of the clustering.

    Description:
        Reports statistics of the clustering - the number and fraction of cells
        in each cluster, gene expression values and dimension reduction plots,
        and, when TCR/BCR data are configured, stats of the TCR clones/clusters
        per cluster. The statistics and plots are produced by `Seurat` and
        `scplotter` (R), and by `clustree` for the `clustrees` plots.

    Base class:
        `biopipen.ns.scrna.SeuratClusterStats`

    Deviations:
        `dimplots` is overridden with the base's
        `{"Dimensional reduction plot": {"label": True}}` entry. When TCR/BCR
        data are configured, a second plot, `VDJ Presence` (grouped by
        `VDJ_Presence`), is added to it; the class body adds that entry only if
        VDJ data are present. `envs_depth` is also set to 3, so that the nested
        `envs` of the plots can be given in the configuration file.
    """

    requires = CombinedInput  # type: ignore
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


@when(VDJInput, requires=CombinedInput)  # type: ignore
class ClonalStats(ClonalStats_):
    """Visualize the clonal information.

    Description:
        Visualizes the clonal information of the TCR/BCR data - clonal volume,
        diversity, overlaps between groups and so on. The plots are drawn by
        `scplotter` (R).

    Base class:
        `biopipen.ns.tcr.ClonalStats`

    Deviations:
        The parameters are passed through to `scplotter` unchanged: this
        process overrides no inherited parameter and adds none.
    """

    envs_depth = 3
    order = 8


@when("ScFGSEA" in config, requires=CombinedInput)
class ScFGSEA(ScFGSEA_):
    """Gene set enrichment analysis for cells in different groups using `fgsea`

    Description:
        Performs gene set enrichment analysis on the expression data for a
        variety of groupings, including ones taken from the metadata and from
        the TCR/BCR data. The testing is done by `fgsea` (R/Bioconductor).

    Base class:
        `biopipen.ns.scrna.ScFGSEA`

    Deviations:
        The parameters are passed through to `fgsea` unchanged: this process
        overrides no inherited parameter and adds none.
    """

    order = 9


@when("PseudoBulkDEG" in config, requires=CombinedInput)
@annotate.format_doc()
class PseudoBulkDEG(PseudoBulkDEG_):
    """{{Summary}}

    SeeAlso:
        - [biopipen.ns.scrna.PseudoBulkDEG](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.PseudoBulkDEG)
        - [ClusterMarkers](./ClusterMarkers.md) for examples of marker and enrichment plots

    Description:
        Performs pseudo-bulk differential gene expression analysis between the
        groups of cells. The testing is done by `DESeq2` (R), the default of
        `envs.tool`, and the results are visualized with `plotthis`/`scplotter`.

    Base class:
        `biopipen.ns.scrna.PseudoBulkDEG`

    Deviations:
        The parameters are passed through unchanged: this process overrides no
        inherited parameter and adds none.
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

    Description:
        Finds the markers between different groups of cells and runs an
        enrichment analysis on them. Markers are found by
        `Seurat::FindMarkers()` and the enrichment is done by `enrichr`.

    Base class:
        `biopipen.ns.scrna.MarkersFinder`

    Deviations:
        The parameters are passed through to `Seurat` unchanged: this process
        overrides no inherited parameter and adds none.
    """  # noqa: E501

    order = 11


@when(VDJInput and "CDR3AAPhyschem" in config, requires=CombinedInput)  # type: ignore
class CDR3AAPhyschem(CDR3AAPhyschem_):
    """CDR3 AA physicochemical feature analysis

    Description:
        Runs a regression between two groups of cells (for example Treg vs
        Tconv) at different lengths of CDR3 amino-acid sequences, for each
        physicochemical feature of the amino acids (hydrophobicity, volume and
        isoelectric point). The modelling is done by `glmnet` (R).

    Base class:
        `biopipen.ns.tcr.CDR3AAPhyschem`

    Deviations:
        The parameters are passed through to `glmnet` unchanged: this process
        overrides no inherited parameter and adds none.
    """

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

    scrna_metabolic_landscape.p_input.requires = CombinedInput  # type: ignore
    scrna_metabolic_landscape.p_input.order = 99
