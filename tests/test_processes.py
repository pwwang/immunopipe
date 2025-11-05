import pytest  # noqa: F401

from .conftest import run_process


@pytest.mark.forked
def test_sampleinfo(request):
    outdir = run_process(
        "SampleInfo", "sampleinfo.config.toml", request=request
    )
    assert outdir.joinpath("Age_distribution-Histogram-.png").is_file()


@pytest.mark.forked
def test_seuratpreparing(request):
    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing.config.toml",
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_seuratpreparing_sct(request):
    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing_sct.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_seuratclustering(request):
    outdir = run_process(
        "SeuratClustering",
        "seuratclustering.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_clustermarkers(request):
    outdir = run_process(
        "ClusterMarkers",
        "clustermarkers.config.toml",
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
        "torbcellselection.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.qs").is_file()


@pytest.mark.forked
def test_screploading(request):
    outdir = run_process(
        "ScRepLoading",
        "screp.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep.qs").is_file()


@pytest.mark.forked
def test_screpcombiningexpression(request):
    outdir = run_process(
        "ScRepCombiningExpression",
        "screp.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep.qs").is_file()


@pytest.mark.forked
def test_seuratclusterstats(request):
    outdir = run_process(
        "SeuratClusterStats",
        "seuratclusterstats.config.toml",
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
        "clonalstats.config.toml",
        export=True,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.scRep.clonalstats").is_dir()
