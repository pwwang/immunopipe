import sys
from pathlib import Path
from argx import Namespace
from simpleconf import Config

from .utils import run_r


def check_dim(args: Namespace) -> None:
    """Check dimensions of the Seurat object before and after QC."""

    signature_file = Path(args.workdir) / "SeuratPreparing" / "0" / "job.signature.toml"
    if not signature_file.exists():
        print(f"Signature file not found: {signature_file}")
        print("Please ensure the pipeline has been run up to SeuratPreparing process.")
        sys.exit(1)

    signature = Config.load_one(signature_file, loader="toml")
    seurat_path = Path(signature["output"]["data"]["outfile"])
    outdir = seurat_path.parent
    check_dim_rfile = Path(__file__).parent / "check_dim.R"
    script = (
        f"source('{check_dim_rfile}')",
        f"check_dim('{outdir}')",
    )
    rc, stdout, stderr = run_r(args.rscript, script)
    if rc != 0:
        print("Error running R script:")
        if stderr:
            print(stderr)
        sys.exit(rc)
    if stdout:
        print(stdout)
