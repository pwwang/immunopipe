import pytest  # noqa: F401

from .conftest import run_process


def test_sampleinfo(tmp_path):

    workdir = run_process("SampleInfo", "sampleinfo.config.toml", tmp_path)
    print(workdir)
