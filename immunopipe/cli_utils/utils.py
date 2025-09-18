from __future__ import annotations
from typing import Sequence


def run_r(
    rscript: str,
    script: str | Sequence[str],
) -> tuple[int, str | None, str | None]:
    """Run an R script using Rscript command.

    Args:
        rscript: Path to Rscript executable.
        script: R script as a string or a list of strings (lines).

    Returns:
        A tuple of (return code, stdout, stderr).
    """
    import subprocess
    import tempfile

    if isinstance(script, str):
        script_content = script
    else:
        script_content = "\n".join(script)

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".R", delete=True) as tf:
        tf.write(script_content)
        tf.flush()
        process = subprocess.run(
            [rscript, tf.name],
            capture_output=True,
            text=True,
        )
        return process.returncode, process.stdout, process.stderr
