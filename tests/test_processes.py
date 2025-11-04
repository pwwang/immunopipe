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
