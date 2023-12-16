"""Process definition"""
from datar.tibble import tibble
from pipen.utils import mark, is_loading_pipeline
from pipen_annotate import annotate
from pipen_args import config
from pipen_filters.filters import FILTERS

# biopipen processes
from biopipen.ns.delim import SampleInfo as SampleInfo_
from biopipen.ns.tcr import (
    ImmunarchLoading as ImmunarchLoading_,
    Immunarch as Immunarch_,
    CloneResidency as CloneResidency_,
    TCRClustering as TCRClustering_,
    TCRClusterStats as TCRClusterStats_,
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
from .inhouse import TCellSelection as TCellSelection_

toml_dumps = FILTERS["toml_dumps"]
just_loading = is_loading_pipeline()
DOC_BASEURL = "https://pwwang.github.io/immunopipe"


@annotate.format_doc(indent=1, vars={"baseurl": DOC_BASEURL})
class SampleInfo(SampleInfo_):
    """{{Summary}}

    This process is the entrance of the pipeline. It just pass by input file and list
    the sample information in the report.

    To specify the input file in the configuration file, use the following

    ```toml
    [SampleInfo.in]
    infile = [ "path/to/sample_info.txt" ]
    ```

    Or with `pipen-board`, find the `SampleInfo` process and click the `Edit` button.
    Then you can specify the input file here

    ![infile]({{baseurl}}/processes/images/SampleInfo-infile.png)

    Theroetically, we can have multiple input files. However, it is not tested yet.
    If you have multiple input files to run, please run it with a different pipeline
    instance (configuration file).

    For the content of the input file, please see details
    [here]({{baseurl}}/preparing-input.md#metadata).

    Once the pipeline is finished, you can see the sample information in the report

    ![report]({{baseurl}}/processes/images/SampleInfo-report.png)

    Note that the required `RNAData` and `TCRData` columns are not shown in the report.
    They are used to specify the paths of the `scRNA-seq` and `scTCR-seq` data, respectively.

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

        | Subject | Sample | Source | Score |
        | ------- | ------ | ------ | ----- |
        | A       | A1     | Tumor  | 1     |
        | A       | A2     | Numor  | 8     |
        | A       | A3     | Tumor  |3      |
        | A       | A4     | Normal |8      |
        | B       | B1     | Tumor  |2      |
        | B       | B2     | Normal |8      |
        | B       | B3     | Tumor  |4      |
        | B       | B4     | Normal |8      |
        | C       | C1     | Tumor  |9      |
        | C       | C2     | Normal |3      |
        | C       | C3     | Tumor  |7      |
        | C       | C4     | Normal |3      |
        | D       | D1     | Tumor  |10     |
        | D       | D2     | Normal |5      |
        | D       | D3     | Tumor  |10     |
        | D       | D4     | Normal |5      |
        | E       | E1     | Tumor  |6      |
        | E       | E2     | Normal |5      |
        | E       | E3     | Tumor  |6      |
        | E       | E4     | Normal |5      |
        | F       | F1     | Tumor  |8      |
        | F       | F2     | Normal |10     |
        | F       | F3     | Tumor  |8      |
        | F       | F4     | Normal |10     |

        ### Count the number of samples per Source

        ```toml
        [SampleInfo.envs.stats]
        Samples_Source = { "group": "Source" }
        Samples_Source_each_Subject = { "group": "Source", "each": "Subject" }
        ```

        ![Samples_Source]({{baseurl}}/processes/images/SampleInfo_Samples_Source.png)

        ### Explore the distribution of the Score

        ```toml
        [SampleInfo.envs.stats.Score_Source_vlnbox]
        on = "Score"
        group = "Source"
        plot = "violin+box"

        [SampleInfo.envs.stats.Score_Source_each_Subject_vlnbox]
        on = "Score"
        group = "Source"
        plot = "violin+box"
        each = "Subject"
        ```

        ![Score_Source]({{baseurl}}/processes/images/SampleInfo_Score_Source.png)

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
    """  # noqa: E501
    envs = {"exclude_cols": "TCRData,RNAData"}


class ImmunarchLoading(ImmunarchLoading_):
    requires = SampleInfo


@annotate.format_doc(indent=1, vars={"baseurl": DOC_BASEURL})
class SeuratPreparing(SeuratPreparing_):
    """{{Summary}}

    See also [Preparing the input](../preparing-input.md#scRNA-seq-data).

    Metadata:
        Here is the demonstration of basic metadata for the `Seurat` object. Future
        processes will use it and/or add more metadata to the `Seurat` object.

        ![SeuratPreparing-metadata]({{baseurl}}/processes/images/SeuratPreparing-metadata.png)
    """  # noqa: E501
    requires = SampleInfo


if "ModuleScoreCalculator" in config or just_loading:
    # Define the process first.
    # We need to setup the connections later
    @annotate.format_doc(indent=2, vars={"baseurl": DOC_BASEURL})
    class ModuleScoreCalculator(ModuleScoreCalculator_):
        """{{Summary}}

        Metadata:
            The metadata of the `Seurat` object will be updated with the module scores:

            ![ModuleScoreCalculator-metadata]({{baseurl}}/processes/images/ModuleScoreCalculator-metadata.png)
        """  # noqa: E501

    if "TCellSelection" not in config and not just_loading:
        ModuleScoreCalculator.requires = SeuratPreparing
        SeuratPreparing = ModuleScoreCalculator


@annotate.format_doc(indent=1)
class SeuratClusteringOfAllCells(SeuratClustering):
    """Cluster all cells using Seurat

    This process will perform clustering on all cells using
    [`Seurat`](https://satijalab.org/seurat/) package.
    The clusters will then be used to select T cells by
    [`TCellSelection`](TCellSelection.md) process.

    {{*Summary.long}}

    /// Note
    If all you cells are T cells, this process will perform clustering on all cells
    and the results will be used on downstream T-cell analyses and the integrative
    analyses. At the same time, you should leave the
    [`TCellSelection`](TCellSelection.md) process out of the pipeline, by not setting
    anything for the process in the configuration file.
    ///
    """
    requires = SeuratPreparing


if "ClusterMarkersOfAllCells" in config or just_loading:
    @annotate.format_doc(indent=2)
    class ClusterMarkersOfAllCells(MarkersFinder_):
        """Markers for clusters of all cells.

        /// Tip | Added in 0.9.0
        `ClusterMarkersOfAllCells` is added in `0.9.0` and is optional by default.
        ///

        See also [ClusterMarkers](./ClusterMarkers.md).

        Envs:
            cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
            each (hidden;readonly): {{Envs.each.help | indent: 12}}.
            ident-1 (hidden;readonly): {{Envs["ident-1"].help | indent: 12}}.
            ident-2 (hidden;readonly): {{Envs["ident-2"].help | indent: 12}}.
            mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
            prefix_each (hidden;readonly): {{Envs.prefix_each.help | indent: 12}}.
            section (hidden;readonly): {{Envs.section.help | indent: 12}}.
        """
        requires = SeuratClusteringOfAllCells
        envs = {
            "cases": {"Cluster": {}},
            "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
        }
        order = 2


if "TopExpressingGenesOfAllCells" in config or just_loading:
    @annotate.format_doc(indent=2)
    class TopExpressingGenesOfAllCells(TopExpressingGenes_):
        """Top expressing genes for clusters of all cells.

        /// Tip | Added in 0.9.0
        `TopExpressingGenesOfAllCells` is added in `0.9.0` and is optional by default.
        ///

        {{*Summary.long}}

        See also [TopExpressingGenes](./TopExpressingGenes.md).

        Envs:
            cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
            each (hidden;readonly): {{Envs.each.help | indent: 12}}.
            group-by (hidden;readonly): {{Envs["group-by"].help | indent: 12}}.
            ident (hidden;readonly): {{Envs.ident.help | indent: 12}}.
            mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
            prefix_each (hidden;readonly): {{Envs.prefix_each.help | indent: 12}}.
            section (hidden;readonly): {{Envs.section.help | indent: 12}}.
        """
        requires = SeuratClusteringOfAllCells
        envs = {"cases": {"Cluster": {}}}
        order = 3


if "TCellSelection" in config or just_loading:
    class TCellSelection(TCellSelection_):
        requires = [SeuratClusteringOfAllCells, ImmunarchLoading]
        input_data = lambda ch1, ch2: tibble(ch1, ch2.iloc[:, 1])

    if "ModuleScoreCalculator" in config or just_loading:
        ModuleScoreCalculator.requires = TCellSelection
        ModuleScoreCalculator.input_data = lambda ch1: ch1.iloc[:, [0]]
        TCellSelection = ModuleScoreCalculator

    @annotate.format_doc(indent=2, vars={"baseurl": DOC_BASEURL})
    class SeuratClusteringOfTCells(SeuratClustering):
        """Cluster the T cells selected by `TCellSelection`.

        If nothing is set for `TCellSelection` in the config file, meaning
        all cells are T cells, this process will be skipped.

        See also: [SeuratClusteringOfAllCells](./SeuratClusteringOfAllCells.md).

        Metadata:
            The metadata of the `Seurat` object will be updated with the cluster
            assignments:

            ![SeuratClusteringOfTCells-metadata]({{baseurl}}/processes/images/SeuratClusteringOfTCells-metadata.png)
        """
        requires = TCellSelection
        input_data = lambda ch1: ch1.iloc[:, [0]]
else:
    SeuratClusteringOfTCells = SeuratClusteringOfAllCells


@annotate.format_doc(indent=1, vars={"baseurl": DOC_BASEURL})
class CellTypeAnnotation(CellTypeAnnotation_):
    """Annotate the T cell clusters.

    {{*Summary}}

    The `<workdir>` is typically `./.pipen` and the `<pipline_name>` is `Immunopipe`
    by default.

    /// Note
    When cell types are annotated, the old `seurat_clusters` column will be renamed
    to `seurat_clusters_id`, and the new `seurat_clusters` column will be added.
    ///

    Metadata:
        When `envs.tool` is `direct` and `envs.cell_types` is empty, the metadata of
        the `Seurat` object will be kept as is.

        When `envs.newcol` is specified, the original `seurat_clusters` column will be
        kept is, and the annotated cell types will be saved in the new column.
        Otherwise, the original `seurat_clusters` column will be replaced by the
        annotated cell types and the original `seurat_clusters` column will be
        saved at `seurat_clusters_id`.

        ![CellTypeAnnotation-metadata]({{baseurl}}/processes/images/CellTypeAnnotation-metadata.png)

    """
    requires = SeuratClusteringOfTCells
    # Change the default to direct, which doesn't do any annotation
    envs = {"tool": "direct", "sctype_db": None}


@annotate.format_doc(indent=1)
class ClusterMarkers(MarkersFinder_):
    """Markers for clusters of T cells.

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

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        group-by (hidden;readonly): {{Envs["group-by"].help | indent: 12}}.
        ident-1 (hidden;readonly): {{Envs["ident-1"].help | indent: 12}}.
        ident-2 (hidden;readonly): {{Envs["ident-2"].help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
        prefix_each (hidden;readonly): {{Envs.prefix_each.help | indent: 12}}.
        section (hidden;readonly): {{Envs.section.help | indent: 12}}.
    """  # noqa: E501
    requires = CellTypeAnnotation
    envs = {
        "cases": {"Cluster": {}},
        "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
    }
    order = 2


@annotate.format_doc(indent=1)
class TopExpressingGenes(TopExpressingGenes_):
    """Top expressing genes for clusters of T cells.

    {{*Summary.long}}

    This process finds the top expressing genes of clusters of T cells, and also
    performs the enrichment analysis against the genes.

    The enrichment analysis is done by [`enrichr`](https://maayanlab.cloud/Enrichr/).

    /// Note
    There are other environment variables also available. However, they should not be
    used in this process. Other environment variables are used for more complicated
    cases for investigating top genes
    (See [`biopipen.ns.scrna.TopExpressingGenes`](https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/#biopipen.ns.scrna.TopExpressingGenes) for more details).

    If you are using `pipen-board` to run the pipeline
    (see [here](../running.md#run-the-pipeline-via-pipen-board) and
    [here](../running.md#run-the-pipeline-via-pipen-board-using-docker-image)),
    you may see the other environment variables of this process are hidden and readonly.
    ///

    Envs:
        cases (hidden;readonly): {{Envs.cases.help | indent: 12}}.
        each (hidden;readonly): {{Envs.each.help | indent: 12}}.
        group-by (hidden;readonly): {{Envs["group-by"].help | indent: 12}}.
        ident (hidden;readonly): {{Envs.ident.help | indent: 12}}.
        mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.
        prefix_each (hidden;readonly): {{Envs.prefix_each.help | indent: 12}}.
        section (hidden;readonly): {{Envs.section.help | indent: 12}}.
    """  # noqa: E501
    requires = CellTypeAnnotation
    envs = {"cases": {"Cluster": {}}}
    order = 3


if "TESSA" in config or just_loading:
    @annotate.format_doc(indent=2, vars={"baseurl": DOC_BASEURL})
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

            ![TESSA-metadata]({{baseurl}}/processes/images/TESSA-metadata.png)
        """
        requires = ImmunarchLoading, CellTypeAnnotation
        input_data = lambda ch1, ch2: tibble(ch1.iloc[:, 1], ch2)

    CellTypeAnnotation = TESSA


@annotate.format_doc(indent=1, vars={"baseurl": DOC_BASEURL})
class SeuratMetadataMutater(SeuratMetadataMutater_):
    """Attach TCR clone information as meta columns to Seurat object

    This process is used to integrate scTCR-seq data into the `Seurat` object.
    The scTCR-seq data is loaded by [ImmunarchLoading](./ImmunarchLoading.md) process.
    The integration is done by matching the barcodes from the `Seurat` object and
    the scTCR-seq data.
    The barcodes from the scTCR-seq data are prefixed with the sample name,
    for example, `Sample1_AAACCTGAGAAGGCTA-1`. The prefix is specified by the `prefix`
    environment variable in the [ImmunarchLoading](./ImmunarchLoading.md) process.

    [ImmunarchLoading](./ImmunarchLoading.md) process will generate a text file with
    the information for each cell.
    `ImmunarchLoading.envs.metacols` can be used to specify the columns to be exported
    to the text file, which will then be integrated into the `Seurat` object
    by this process.

    You may also use `envs.mutaters` to add new columns to the metadata.
    These columns can be used for downstream analysis.
    An additional column `TCR_Presence` is added so later on we can overlay the
    TCR presence on the UMAP plot in
    [`SeuratClusteringOfTCells`](./SeuratClusteringOfTCells.md) process.

    /// Warning
    If you are modifying `envs.mutaters`, make sure you keep the `TCR_Presence` column.
    Because by default, `SeuratClusteringOfTCells` process will use this column to
    overlay the TCR presence on the UMAP plot.
    ///

    {{*Summary.long}}

    Metadata:
        The metadata of the `Seurat` object will be updated with information from TCR
        data:

        ![SeuratMetadataMutater-metadata]({{baseurl}}/processes/images/SeuratMetadataMutater-metadata.png)

        All of the columns above can be used for downstream analysis.
    """
    requires = CellTypeAnnotation, ImmunarchLoading
    input_data = lambda ch1, ch2: tibble(ch1.iloc[:, 0], ch2.iloc[:, 1])
    envs = {
        "mutaters": {
            "TCR_Presence": 'if_else(is.na(CDR3.aa), "TCR_absent", "TCR_present")'
        }
    }


if (
    "TCRClustering" in config
    or "TCRClusterStats" in config
    or just_loading
):
    @annotate.format_doc(indent=2)
    class TCRClustering(TCRClustering_):
        """{{Summary.short}}

        You can disable this by remving the whole sections of
        TCRClustering and TCRClusterStats in the config file.

        {{*Summary.long}}
        """
        requires = ImmunarchLoading
        input_data = lambda ch1: ch1.iloc[:, [0]]

    @mark(board_config_hidden=True)
    @annotate.format_doc(indent=2, vars={"baseurl": DOC_BASEURL})
    class TCRClusters2Seurat(SeuratMetadataMutater_):
        """Attach TCR clusters as meta columns to Seurat object

        {{*Summary.long}}

        This process is used to merge the cluster assignments from
        [TCRClustering](./TCRClustering.md) to the `Seurat` object.
        The cluster assignments are prefixed with `S_` or `M_` to indicate whether
        a cluster has only one unique CDR3 sequence or multiple CDR3 sequences.
        Note that a cluster with `S_` prefix may still have multiple cells,
        as the same CDR3 sequence may be shared by multiple cells.
        The cluster assignments are saved in the `Seurat` object at `TCR_Cluster`
        column in `seurat_object@meta.data` in `R`.

        Other two columns are also added to the `Seurat` object: `TCR_Cluster_Size`
        and `TCR_Cluster_Size1`. The `TCR_Cluster_Size` column contains the number of
        cells in each cluster, while the `TCR_Cluster_Size1` column contains the
        number of unique CDR3 sequences in each cluster.

        Those columns can be then used for further downstream analysis. For example,
        you can find the markers for the TCR cluster (i.e. `S_1` vs `S_2`) in
        each seurat cluster by

        ```toml
        [MarkersFinder.envs]
        group-by = "TCR_Cluster"
        ident-1 = "S_1"
        ident-2 = "S_2"
        each = "seurat_clusters"
        ```

        Envs:
            mutaters (hidden;readonly): {{Envs.mutaters.help | indent: 12}}.

        Metadata:
            The metadata of the `Seurat` object will be updated with the TCR cluster
            assignments and their sizes:

            ![TCRClusters2Seurat-metadata]({{baseurl}}/processes/images/TCRClusters2Seurat-metadata.png)
        """
        requires = SeuratMetadataMutater, TCRClustering
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1.rdsfile, metafile=ch2.clusterfile
        )

    class TCRClusterStats(TCRClusterStats_):
        requires = TCRClustering
        input_data = lambda ch1: ch1.iloc[:, [0]]

else:
    TCRClusters2Seurat = SeuratMetadataMutater


class SeuratClusterStats(SeuratClusterStats_):
    requires = TCRClusters2Seurat
    order = -1
    envs = {
        "dimplots": {
            "Dimensional reduction plot": {
                "label": True,
                "label-box": True,
                "repel": True,
            },
            "TCR presence": {
                "ident": "TCR_Presence",
                "order": "TCR_absent",
                "cols": ["#FF000066", "gray"],
            },
        },
    }


@annotate.format_doc(indent=1)
class Immunarch(Immunarch_):
    """{{Summary.short}}

    /// Tip | Changed in `0.10.0`
    `envs.mutaters` are now applied at cell level.

    Seurat clustering information and other information are added at cell level, which
    can be used to subset the cells for listed analyses.

    You can now use `subset` to subset the cells for listed analyses, at cell level.
    ///

    {{*Summary.long | str |
    replace: '!!#biopipennstcrimmunarchloading', './ImmunarchLoading.md'}}
    """
    requires = ImmunarchLoading, TCRClusters2Seurat
    input_data = lambda ch1, ch2: tibble(
        immdata=ch1.iloc[:, 0],
        metafile=ch2.iloc[:, 0],
    )


if "CellsDistribution" in config or just_loading:
    class CellsDistribution(CellsDistribution_):
        requires = TCRClusters2Seurat
        order = 8


if "CloneResidency" in config or just_loading:
    class CloneResidency(CloneResidency_):
        requires = ImmunarchLoading, TCRClusters2Seurat
        input_data = lambda ch1, ch2: tibble(
            immdata=ch1.iloc[:, 0],
            metafile=ch2.iloc[:, 0],
        )
        order = -3


# if "CloneHeterogeneity" in config or just_loading:
#     class CloneHeterogeneity(CloneHeterogeneity):
#         requires = SeuratMetadataMutater


if "RadarPlots" in config or just_loading:
    @annotate.format_doc(indent=2)
    class RadarPlots(RadarPlots_):
        """{{Summary}}

        Envs:
            mutaters: {{Envs.mutaters.help | indent: 12}}.
                See also
                [`mutating the metadata`](../configurations.md#mutating-the-metadata).
        """
        requires = TCRClusters2Seurat
        order = -2


if "ScFGSEA" in config or just_loading:
    class ScFGSEA(ScFGSEA_):
        requires = TCRClusters2Seurat
        order = 4


if "MarkersFinder" in config or just_loading:
    @annotate.format_doc(indent=2)
    class MarkersFinder(MarkersFinder_):
        """{{Summary.short}}

        `MarkersFinder` is a process that wraps the
        [`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers)
        function, and performs enrichment analysis for the markers found.

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

            By default, `group-by` is `seurat_clusters`, and `ident-1` and `ident-2`
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

            `each` is used to separate the cells into different cases. `group-by`
            is still `seurat_clusters`.

            ```toml
            [<Proc>.envs]
            group-by = "seurat_clusters"
            each = "Group"
            ```

            - A:Cluster
                - 1 (vs 2)
                - 2 (vs 1)
            - B:Cluster
                - 3 (vs 4)
                - 4 (vs 3)

            ### With `ident-1` only

            `ident-1` is used to specify the first group of cells to compare.
            Then the rest of the cells in the case are used for `ident-2`.

            ```toml
            [<Proc>.envs]
            group-by = "seurat_clusters"
            ident-1 = "1"
            ```

            - Cluster
                - 1 (vs 2, 3, 4)

            ### With both `ident-1` and `ident-2`

            `ident-1` and `ident-2` are used to specify the two groups of cells to
            compare.

            ```toml
            [<Proc>.envs]
            group-by = "seurat_clusters"
            ident-1 = "1"
            ident-2 = "2"
            ```

            - Cluster
                - 1 (vs 2)

            ### Multiple cases

            ```toml
            [<Proc>.envs.cases]
            c1_vs_c2 = {ident-1 = "1", ident-2 = "2"}
            c3_vs_c4 = {ident-1 = "3", ident-2 = "4"}
            ```

            - DEFAULT:c1_vs_c2
                - 1 (vs 2)
            - DEFAULT:c3_vs_c4
                - 3 (vs 4)

            The `DEFAULT` section name will be ignored in the report. You can specify
            a section name other than `DEFAULT` for each case to group them
            in the report.
        """
        requires = TCRClusters2Seurat
        envs = {"assay": "RNA"}
        order = 5


if "MetaMarkers" in config or just_loading:
    class MetaMarkers(MetaMarkers_):
        requires = TCRClusters2Seurat
        order = 6


if "CDR3AAPhyschem" in config or just_loading:
    class CDR3AAPhyschem(CDR3AAPhyschem_):
        requires = ImmunarchLoading, TCRClusters2Seurat
        input_data = lambda ch1, ch2: tibble(
            immdata=ch1.rdsfile,
            srtobj=ch2.rdsfile,
        )
        order = 9


if "ScrnaMetabolicLandscape" in config or just_loading:
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
    scrna_metabolic_landscape.p_features_intra_subset.__doc__ = (
        scrna_metabolic_landscape.p_features_intra_subset.__doc__.replace(
            "!!#biopipennsscrna_metabolic_landscapemetabolicfeatures",
            "./MetabolicFeatures.md",
        )
    )
