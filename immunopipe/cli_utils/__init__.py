from argx import ArgumentParser


def main(argv: list[str]) -> None:
    """Entry point for the CLI utility functions."""
    parser = ArgumentParser(
        prog="immunopipe utils",
        description="Utility commands to assist with immunopipe tasks"
    )
    gene_command = parser.add_command(
        "check-genes",
        help=(
            "Check exsistence of gene symbols in the data.\n"
            "This is useful for verifying if the gene symbols in your dataset "
            "when visualizing their expressions in SeuratClusterStats\n"
        )
    )
    gene_command.add_argument(
        "-w",
        "--workdir",
        required=True,
        type=str,
        help="Working directory of the pipeline. Typically ./pipen/<pipeline name>",
    )
    gene_command.add_argument(
        "--rscript",
        type=str,
        default="Rscript",
        help="Path to Rscript executable to run the R script",
    )
    gene_command.add_argument(
        "-g",
        "--genes",
        type=str,
        required=True,
        help=(
            "Comma-separated list of gene symbols to check or "
            "a file path starting with 'file://' containing gene symbols (one per line)"
        )
    )
    gene_command.add_argument(
        "--assay",
        type=str,
        help=(
            "Assay name to check the genes against. "
            "The seurat object will be pulled from SeuratPreparing process."
        )
    )

    dim_command = parser.add_command(
        "check-dim",
        help=(
            "Check dimensions of the Seurat object before and after QC.\n"
            "This is useful for verifying the effect of QC filtering steps.\n"
        )
    )
    dim_command.add_argument(
        "-w",
        "--workdir",
        required=True,
        type=str,
        help="Working directory of the pipeline. Typically ./pipen/<pipeline name>",
    )
    dim_command.add_argument(
        "--rscript",
        type=str,
        default="Rscript",
        help="Path to Rscript executable to run the R script",
    )

    args = parser.parse_args(argv)
    if args.COMMAND == "check-genes":
        from .check_genes import check_genes
        check_genes(args)
    elif args.COMMAND == "check-dim":
        from .check_dim import check_dim
        check_dim(args)
