import pytest  # noqa: F401

from .conftest import run_process

pytest_order = 2


@pytest.mark.forked
def test_immunarchloading(tmp_path, request):

    outdir = run_process(
        "ImmunarchLoading",
        "immunarchloading.config.toml",
        tmp_path,
        request=request,
        export=False,
    )

    assert outdir.joinpath("sampleinfo.immunarch.RDS").is_file()
    assert outdir.joinpath("sampleinfo.tcr.txt").is_file()
