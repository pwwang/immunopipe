import pytest  # noqa: F401

from .conftest import run_process

pytest_order = 1


@pytest.mark.forked
def test_sampleinfo(tmp_path, request):
    outdir = run_process(
        "SampleInfo", "sampleinfo.config.toml", tmp_path, request=request
    )
    assert outdir.joinpath("Age_distribution-Histogram-.png").is_file()
