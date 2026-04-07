"""
Process docker/cleanup.list to remove build-time artefacts from /opt/conda.

Manifest format (see cleanup.list for annotated examples):
  path/to/remove     — exact path relative to /opt/conda; removed silently if present
  path/*/glob        — shell glob (* and ** supported); all matches are removed
  !path/to/keep      — exclusion: never remove this path even if matched above
  # comment          — ignored
  (blank line)       — ignored
"""

import glob as _glob
import os
import shutil

BASE = "/opt/conda"
HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "cleanup.list")


def _load_manifest():
    exclusions = []
    patterns = []
    with open(MANIFEST) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("!"):
                exclusions.append(os.path.join(BASE, line[1:]))
            else:
                patterns.append(line)
    return exclusions, patterns


def _is_excluded(path, exclusions):
    for exc in exclusions:
        if path == exc or path.startswith(exc + os.sep):
            return True
    return False


def _remove(path, exclusions):
    if _is_excluded(path, exclusions):
        return
    try:
        if os.path.islink(path) or os.path.isfile(path):
            os.unlink(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
    except OSError:
        pass


def main():
    exclusions, patterns = _load_manifest()
    for pattern in patterns:
        full = os.path.join(BASE, pattern)
        if "*" in pattern or "?" in pattern:
            # Sort by length descending so children are removed before parents,
            # preventing "already removed" errors when both parent and child match.
            matches = sorted(
                _glob.glob(full, recursive=True), key=len, reverse=True
            )
            for p in matches:
                _remove(p, exclusions)
        else:
            _remove(full, exclusions)


if __name__ == "__main__":
    main()
