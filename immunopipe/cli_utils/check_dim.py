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
    if (
        not seurat_path.exists()
        and seurat_path.as_posix().startswith("/mnt/disks/.cwd")
    ):
        # The pipeline is likely run in a container or google batch jobs
        # Conver the path to a local path based on the path to the signature file
        signature_file = signature_file.resolve()
        cwd = signature_file.parents[4]
        seurat_path = cwd / seurat_path.relative_to("/mnt/disks/.cwd")
        if not seurat_path.exists() and seurat_path.parent.is_symlink():
            seurat_path = seurat_path.parent.resolve() / seurat_path.name

        if seurat_path.as_posix().startswith("/mnt/disks/.cwd"):
            seurat_path = cwd / seurat_path.relative_to("/mnt/disks/.cwd")

    if not seurat_path.exists():
        print(f"Seurat file not found: {seurat_path}")
        sys.exit(1)

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
