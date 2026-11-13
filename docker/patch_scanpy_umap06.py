#!/usr/bin/env python
"""Restore scanpy compatibility with umap-learn >= 0.6.

Backport of the two code changes of scverse/scanpy#4480 ("fix: support
umap-learn 0.6") onto the scanpy installed in this image:

* ``scanpy/neighbors/_connectivity.py`` -- pass C-contiguous ``knn_indices`` /
  ``knn_dists`` to ``umap.fuzzy_simplicial_set``.  umap 0.6.0 stopped copying
  them internally (umap PR #1275), so scanpy's ``distances[:, :n_neighbors]``
  slice hit numba with a non-contiguous array and raised
  ``No conversion from array(float32, 2d, A) to array(float32, 2d, C) for
  'knn_dists'`` (lmcinnes/umap#1335).  This is the crash that killed
  ``CellTypeAnnotation`` (celltypist over-clustering) on >= 8192 cells.
* ``scanpy/tools/_ingest.py`` -- ask ``UMAP`` for the legacy (0.5.x) layout
  via ``compatibility_layout=True`` when the installed umap knows that
  keyword.

Only the parts that matter for running the pipeline are ported: the upstream
test-file hunks are test-only and are not used here.

Idempotent: a file that already carries the fix is left untouched, so this
script can stay in the build until scanpy >= 1.12.6 is in the image, and then
simply be dropped.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

CONNECTIVITY = Path("neighbors/_connectivity.py")
INGEST = Path("tools/_ingest.py")

MARK_CONTIG = "ascontiguousarray(knn_indices)"
MARK_COMPAT = "compatibility_layout"

_DISTS = re.compile(
    r"(?m)^(?P<i>[ \t]*)knn_indices=knn_indices,[ \t]*\r?\n"
    r"(?P=i)knn_dists=knn_dists,[ \t]*\r?\n"
)


def patch_connectivity(path: Path) -> str:
    src = path.read_text()
    if MARK_CONTIG in src:
        return "already patched"
    if "import numpy as np" not in src:
        raise SystemExit(f"ERROR: {path}: numpy is not imported; refusing to patch")
    m = _DISTS.search(src)
    if m is None:
        raise SystemExit(
            f"ERROR: {path}: could not find the knn_indices/knn_dists kwargs of "
            "umap.fuzzy_simplicial_set; scanpy's code moved -- re-check this patch"
        )
    indent = m.group("i")
    repl = (
        f"{indent}# umap >= 0.6 requires C-contiguous arrays\n"
        f"{indent}knn_indices=np.ascontiguousarray(knn_indices),\n"
        f"{indent}knn_dists=np.ascontiguousarray(knn_dists),\n"
    )
    path.write_text(src[: m.start()] + repl + src[m.end() :])
    return "patched"


def patch_ingest(path: Path) -> str:
    src = path.read_text()
    if MARK_COMPAT in src:
        return "already patched"

    if "from inspect import signature" not in src:
        anchor = "from collections.abc import MutableMapping\n"
        if anchor not in src:
            raise SystemExit(
                f"ERROR: {path}: no anchor to add `from inspect import signature`"
            )
        src = src.replace(anchor, anchor + "from inspect import signature\n", 1)

    lines = src.splitlines(keepends=True)
    start = next(
        (i for i, ln in enumerate(lines) if ln.lstrip().startswith("def _init_umap(")),
        None,
    )
    if start is None:
        raise SystemExit(f"ERROR: {path}: `def _init_umap(` not found")
    call = next(
        (
            i
            for i in range(start, len(lines))
            if lines[i].lstrip().startswith("self._umap = UMAP(")
        ),
        None,
    )
    if call is None:
        raise SystemExit(f"ERROR: {path}: `self._umap = UMAP(` not found in _init_umap")
    if call + 1 >= len(lines) or not lines[call + 1].lstrip():
        raise SystemExit(
            f"ERROR: {path}: `self._umap = UMAP(` is not followed by an indented "
            "argument line; scanpy's code moved -- re-check this patch"
        )

    call_indent = lines[call][: len(lines[call]) - len(lines[call].lstrip())]
    arg_indent = lines[call + 1][: len(lines[call + 1]) - len(lines[call + 1].lstrip())]
    compat = (
        f"{call_indent}# we replicate the legacy (0.5.x) transform pipeline below\n"
        f"{call_indent}compat = (\n"
        f'{call_indent}    dict(compatibility_layout=True)\n'
        f'{call_indent}    if "compatibility_layout" in signature(UMAP).parameters\n'
        f"{call_indent}    else {{}}\n"
        f"{call_indent})\n"
    )
    lines[call:call] = [compat]
    lines.insert(call + 2, f"{arg_indent}**compat,\n")
    path.write_text("".join(lines))
    return "patched"


def main(argv: list[str]) -> int:
    if len(argv) > 1:
        root = Path(argv[1])
    else:
        import scanpy

        root = Path(scanpy.__file__).parent

    if not (root / CONNECTIVITY).is_file():
        raise SystemExit(f"ERROR: {root / CONNECTIVITY} does not exist")
    print(f"scanpy at {root}")
    applied = False
    for rel, fn in ((CONNECTIVITY, patch_connectivity), (INGEST, patch_ingest)):
        state = fn(root / rel)
        applied = applied or state == "patched"
        print(f"  {rel}: {state}")

    print("RESULT:", "patch applied" if applied else "nothing to do (already fixed)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
