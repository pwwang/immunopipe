import pytest  # noqa: F401

from .conftest import run_process


def test_sampleinfo(tmp_path, request):

    outdir = run_process(
        "SampleInfo",
        "sampleinfo.config.toml",
        tmp_path,
        request=request
    )
    assert outdir.joinpath("Age_distribution (boxplot).png").is_file()
    assert outdir.joinpath('Age_distribution_in_each_Diagnosis (violin + boxplot).png').is_file()
    assert outdir.joinpath('Age_distribution_per_Diagnosis (violin + boxplot).png').is_file()
    assert outdir.joinpath('Age_distribution_per_Sex_in_each_Diagnosis (boxplot).png').is_file()
    assert outdir.joinpath('N_Samples_per_Diagnosis (bar).png').is_file()
    assert outdir.joinpath('N_Samples_per_Diagnosis (pie).png').is_file()
