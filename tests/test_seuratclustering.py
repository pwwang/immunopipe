import pytest  # noqa: F401

from .conftest import run_process

pytest_order = 3


@pytest.mark.forked
def test_seuratclusteringofallcells(tmp_path, request):

    outdir = run_process(
        "SeuratClusteringOfAllCells",
        "seuratclusteringofallcells.config.toml",
        tmp_path,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()


@pytest.mark.forked
def test_seuratclustering(tmp_path, request):

    outdir = run_process(
        "SeuratClustering",
        "seuratclustering.config.toml",
        tmp_path,
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()
