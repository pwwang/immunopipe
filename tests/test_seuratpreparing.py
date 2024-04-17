import pytest  # noqa: F401

from .conftest import run_process

pytest_order = 2


@pytest.mark.forked
def test_seuratpreparing(tmp_path, request):

    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing.config.toml",
        tmp_path,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()
    assert outdir.joinpath("DoubletFinder_doublets_singlets.txt").is_file()
    assert outdir.joinpath("DoubletFinder_summary.txt").is_file()
    assert outdir.joinpath("plots", "dim.txt").is_file()
    assert outdir.joinpath("plots", "nCount_RNA.vln.png").is_file()
    assert outdir.joinpath("plots", "nFeature_RNA-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "nFeature_RNA.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-hb-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-hb.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-mt-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-mt.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-plat-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-plat.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-ribo-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-ribo.vln.png").is_file()


@pytest.mark.forked
def test_seuratpreparing_sct(tmp_path, request):

    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing_sct.config.toml",
        tmp_path,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()
    assert outdir.joinpath("plots", "dim.txt").is_file()
    assert outdir.joinpath("plots", "nCount_RNA.vln.png").is_file()
    assert outdir.joinpath("plots", "nFeature_RNA-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "nFeature_RNA.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-hb-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-hb.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-mt-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-mt.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-plat-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-plat.vln.png").is_file()
    assert outdir.joinpath("plots", "percent-ribo-nCount_RNA.scatter.png").is_file()
    assert outdir.joinpath("plots", "percent-ribo.vln.png").is_file()
