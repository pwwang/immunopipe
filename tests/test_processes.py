import pytest  # noqa: F401

from .conftest import run_process


@pytest.mark.forked
def test_sampleinfo(request):
    outdir = run_process(
        "SampleInfo", "SampleInfo.config.toml", request=request
    )
    assert outdir.joinpath("Age_distribution-Histogram-.png").is_file()


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
        "sampleinfo.markers/Cluster/seurat_clusters-All-Markers-/"
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
