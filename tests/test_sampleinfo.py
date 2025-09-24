import pytest  # noqa: F401

from .conftest import run_process

pytestmark = pytest.mark.order(1)


@pytest.mark.forked
def test_sampleinfo(request):
    outdir = run_process(
        "SampleInfo", "sampleinfo.config.toml", request=request
    )
    assert outdir.joinpath("Age_distribution-Histogram-.png").is_file()
