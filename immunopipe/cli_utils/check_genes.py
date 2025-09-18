import sys
from pathlib import Path
from argx import Namespace
from simpleconf import Config

from .utils import run_r


def check_genes(args: Namespace) -> None:
    """Check existence of gene symbols in the data."""

    genes_input = args.genes
    if genes_input.startswith("file://"):
        file_path = genes_input[len("file://") :]
        try:
            with open(file_path, "r") as f:
                genes = [line.strip().split("\t")[0] for line in f if line.strip()]
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            sys.exit(1)
    else:
        genes = [gene.strip() for gene in genes_input.split(",") if gene.strip()]

    genes = list(set(genes))  # Remove duplicates
    genes = [repr(gene) for gene in genes if gene]  # Remove empty strings
    genes = ", ".join(genes)
    genes = f"c({genes})"

    signature_file = Path(args.workdir) / "SeuratPreparing" / "0" / "job.signature.toml"
    if not signature_file.exists():
        print(f"Signature file not found: {signature_file}")
        print("Please ensure the pipeline has been run up to SeuratPreparing process.")
        sys.exit(1)

    signature = Config.load_one(signature_file, loader="toml")
    seurat_path = signature["output"]["data"]["outfile"]
    assay = args.assay
    check_genes_rfile = Path(__file__).parent / "check_genes.R"
    script = (
        f"source('{check_genes_rfile}')",
        f"genes <- {genes}",
        f"assay <- {repr(assay) if assay else 'NULL'}",
        f"check_genes('{seurat_path}', genes=genes, assay=assay)",
    )
    rc, stdout, stderr = run_r(args.rscript, script)
    if rc != 0:
        print("Error running R script:")
        if stderr:
            print(stderr)
        sys.exit(rc)
    if stdout:
        print(stdout)
