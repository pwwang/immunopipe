import pytest  # noqa: F401

from .conftest import CONFIGDIR, run_process, with_argv


def test_seuratmap2ref_is_marked_deprecated():
    from pipen.utils import LOADING_ARGV0

    with with_argv(
        [
            LOADING_ARGV0,
            f"@{CONFIGDIR / 'celltypeannotation_seuratmap2ref.config.toml'}",
        ]
    ):
        from importlib import reload
        from immunopipe import processes

        reload(processes)
        SeuratMap2Ref = processes.SeuratMap2Ref

    assert SeuratMap2Ref.__meta__["deprecated"] == (
        "Use CellTypeAnnotation with envs.tool = 'mapquery' instead."
    )


def test_seuratmap2ref_can_be_disabled(tmp_path):
    configfile = tmp_path / "config.toml"
    configfile.write_text("[SampleInfo]\n")
    with with_argv(["immunopipe", f"@{configfile}"]):
        from importlib import reload
        from immunopipe import processes

        reload(processes)
        assert processes.SeuratMap2Ref is None


@pytest.mark.forked
def test_sampleinfo(request):
    outdir = run_process(
        "SampleInfo", "SampleInfo.config.toml", request=request
    )
    assert outdir.joinpath("Age_distribution-Histogram.png").is_file()


@pytest.mark.forked
def test_seuratpreparing(request):
    outdir = run_process(
        "SeuratPreparing",
        "SeuratPreparing.config.toml",
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_seuratpreparing_sct(request):
    outdir = run_process(
        "SeuratPreparing",
        "SeuratPreparing_sct.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_seuratclustering(request):
    outdir = run_process(
        "SeuratClustering",
        "SeuratClustering.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_clustermarkers(request):
    outdir = run_process(
        "ClusterMarkers",
        "ClusterMarkers.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.markers/Cluster/seurat_clusters-All-Markers/"
        "Top-10-markers-of-all-clusters.png"
    ).is_file()


@pytest.mark.forked
def test_torbcellselection(request):
    outdir = run_process(
        "TOrBCellSelection",
        "TOrBCellSelection.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_torbcellselection_with_cellbarcodes(request):
    outdir = run_process(
        "TOrBCellSelection",
        "TOrBCellSelectionWithCellbarcodes.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_screploading(request):
    outdir = run_process(
        "ScRepLoading",
        "ScRep.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep.qs").is_file()


@pytest.mark.forked
def test_screpcombiningexpression(request):
    outdir = run_process(
        "ScRepCombiningExpression",
        "ScRep.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep.qs").is_file()


@pytest.mark.forked
def test_seuratclusterstats(request):
    outdir = run_process(
        "SeuratClusterStats",
        "SeuratClusterStats.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.scRep.cluster_stats/clustrees/seurat_clusters.clustree.png"
    ).is_file()


@pytest.mark.forked
def test_clonalstats(request):
    outdir = run_process(
        "ClonalStats",
        "ClonalStats.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep.clonalstats").is_dir()
    assert outdir.joinpath(
        "sampleinfo.scRep.clonalstats/Clonal-Statistics/"
        "TRB-Gene-Call-Combinations-Across-Diagnoses.png"
    ).is_file()


@pytest.mark.forked
def test_scfgsea(request):
    outdir = run_process(
        "ScFGSEA",
        "ScFGSEA.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.fgsea").is_dir()


@pytest.mark.forked
def test_ccc(request):
    outdir = run_process(
        "CellCellCommunication",
        "CellCellCommunication.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep-ccc.txt").is_file()


@pytest.mark.forked
def test_cccplots(request):
    outdir = run_process(
        "CellCellCommunicationPlots",
        "CellCellCommunication.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.scRep-ccc_plots/Cell-Cell-Communication-Circos-Plot.png"
    ).is_file()


@pytest.mark.forked
def test_topexpressinggenes(request):
    outdir = run_process(
        "TopExpressingGenes",
        "TopExpressingGenes.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.scRep.top_expressing_genes/Cluster/c1/"
        "enrich.MSigDB_Hallmark_2020.Bar-Plot.png"
    ).is_file()


@pytest.mark.forked
def test_pseudobulkdeg(request):
    outdir = run_process(
        "PseudoBulkDEG",
        "PseudoBulkDEG.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.scRep.pseudobulk_deg/DEG-Analysis/degs.Volcano-Plot.png"
    ).is_file()


@pytest.mark.forked
def test_pseudobulkdeg_edger(request):
    outdir = run_process(
        "PseudoBulkDEG",
        "PseudoBulkDEG_edgeR.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.scRep.pseudobulk_deg/DEG-Analysis/degs.Volcano-Plot.png"
    ).is_file()


@pytest.mark.forked
def test_cdr3aaphyschem(request):
    outdir = run_process(
        "CDR3AAPhyschem",
        "CDR3AAPhyschem.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath(
        "sampleinfo.scRep.cdr3aaphyschem/ALL/estimated_coefficients.png"
    ).is_file()
