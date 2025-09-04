from __future__ import annotations
from typing import Generator

import sys
from contextlib import contextmanager
from pathlib import Path

import pytest
from pipen import Pipen, Proc

RUNNINGDIR = Path(__file__).parent / "running"

DATADIR = RUNNINGDIR / "data"
CONFIGDIR = RUNNINGDIR / "configs"
OUTDIR = RUNNINGDIR / "output"
INTERMEDIATEDIR = RUNNINGDIR / "intermediate"
WORKDIR = RUNNINGDIR / "workdir"


@contextmanager
def with_argv(argv: list[str]) -> Generator:
    """Temporarily change sys.argv."""
    orig_argv = sys.argv
    sys.argv = argv
    yield
    sys.argv = orig_argv


def run_process(
    process: str,
    configfile: str,
    tmp_path: Path = WORKDIR,
    request: pytest.Fixture = None,
    export: bool = True,
    **kwargs,
) -> Path:
    """Run a process with a given config file.

    Args:
        process: The process to test
        configfile: The config file to use
        tmp_path: The temporary directory to use
        request: The pytest request object
        **kwargs: The arguments to pass to the process

    Returns:
        The working directory of the process
    """
    configfile = str(CONFIGDIR / configfile)

    with with_argv(["@pipen"]):
        from immunopipe.pipeline import Immunopipe

    full_pipe = Immunopipe()
    full_pipe.build_proc_relationships()

    for proc in full_pipe.procs:
        if proc.name == process:
            break
    else:
        raise ValueError(f"Process {process} not found in the pipeline.")

    # detech dependent procs
    proc = Proc.from_proc(proc, name=process, **kwargs)

    class Pipeline(Pipen):
        starts = proc
        plugin_opts = {"args_flatten": False}
        if request:
            name = request.node.name[5:]
        else:
            name = process.lower()

        outdir = (OUTDIR if export else INTERMEDIATEDIR) / name

    with with_argv(["immunopipe", f"@{configfile}"]):
        pipe = Pipeline(
            plugins=["-report", "-diagram"],
            workdir=tmp_path,
        )
        if not pipe.run():
            raise RuntimeError("Failed to run the process.")

    return pipe.outdir / process
