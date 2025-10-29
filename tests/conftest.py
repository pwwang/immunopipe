from __future__ import annotations
from typing import Generator

import re
import sys
from contextlib import contextmanager
from pathlib import Path
from importlib import reload
from subprocess import Popen, PIPE

import pytest
from pipen import Pipen, Proc
from pipen.utils import LOADING_ARGV0

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

    # with with_argv(["@pipen"]):
    #     from immunopipe.pipeline import Immunopipe

    # full_pipe = Immunopipe()
    # full_pipe.build_proc_relationships()
    # print(full_pipe.procs)

    # for proc in full_pipe.procs:
    #     if proc.name == process:
    #         break
    # else:
    #     raise ValueError(f"Process {process} not found in the pipeline.")
    with with_argv([LOADING_ARGV0, f"@{configfile}"]):
        # This is loaded already by pytest
        # In order to load all process, we need to reload processes module
        from immunopipe import processes
        from immunopipe import pipeline

        reload(processes)
        reload(pipeline)
        from immunopipe.pipeline import Immunopipe

        pipe = Immunopipe()
        pipe.build_proc_relationships()
        proc = [p for p in pipe.procs if p.name == process][0]
        proc.nexts = []

    # detech dependent procs
    proc = Proc.from_proc(proc, name=process, **kwargs)
    Pipen.SETUP = True

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


def dry_run(*args) -> str:
    """Run a dry run with given arguments.

    Args:
        *args: The arguments to pass to the dry run

    Returns:
        The stdout (logs) of the dry run
    """
    command = [
        sys.executable,
        "-m",
        "immunopipe",
        *args,
        "--plugins=+dry",
        "--plugins=-report",
        "--plugins=-log2file",
        "--plugins=-board",
        "--plugins=-poplog",
        "--plugins=-verbose",
        "--plugins=-deprecated",
        "--scheduler=dry",
    ]
    process = Popen(
        command,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
        bufsize=1,
        encoding="utf-8",
    )
    process.wait()
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(
            "".join(
                [
                    "Dry run failed.\n\n",
                    "Command: \n",
                    " ".join(command),
                    "\n\n",
                    "Stdout:\n",
                    stdout,
                    "\n\n",
                    "Stderr:\n",
                    stderr,
                    "\n",
                ]
            )
        )

    stdout = re.sub(r"[\n\s]+", " ", stdout)
    return stdout
