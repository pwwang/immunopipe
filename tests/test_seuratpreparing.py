import pytest  # noqa: F401

from .conftest import run_process

pytest_order = 2


@pytest.mark.forked
def test_seuratpreparing(request):

    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing.config.toml",
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()


@pytest.mark.forked
def test_seuratpreparing_sct(request):

    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing_sct.config.toml",
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()
