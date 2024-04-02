from __future__ import annotations
from typing import Generator

import sys
from contextlib import contextmanager
from pathlib import Path

from pipen import Pipen, Proc

DATADIR = Path(__file__).parent / "data"
CONFIGDIR = Path(__file__).parent / "configs"


@contextmanager
def with_argv(argv: list[str]) -> Generator:
    """Temporarily change sys.argv."""
    orig_argv = sys.argv
    sys.argv = argv
    yield
    sys.argv = orig_argv


def run_process(process: str, configfile: str, tmp_path: Path) -> Path:
    """Run a process with a given config file.

    Args:
        process: The process to test
        configfile: The config file to use
        tmp_path: The temporary directory to use

    Returns:
        The working directory of the process
    """
    configfile = str(CONFIGDIR / configfile)

    with with_argv(["@pipen"]):
        from immunopipe import processes

    with with_argv(["immunopipe", f"@{configfile}"]):
        proc = getattr(processes, process)
        # detech dependent procs
        proc = Proc.from_proc(proc, name=process)

        class Pipeline(Pipen):
            starts = proc
            workdir = tmp_path
            plugin_opts = {"args_flatten": False}

        pipe = Pipeline(plugins=["-report", "-diagram"])
        pipe.run()

    return pipe.workdir
