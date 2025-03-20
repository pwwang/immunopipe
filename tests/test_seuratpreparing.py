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


@pytest.mark.forked
def test_seuratpreparing_sct(tmp_path, request):

    outdir = run_process(
        "SeuratPreparing",
        "seuratpreparing_sct.config.toml",
        tmp_path,
        export=False,
        request=request,
    )
    assert outdir.joinpath("sampleinfo.seurat.RDS").is_file()
