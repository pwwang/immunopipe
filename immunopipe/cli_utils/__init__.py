from argx import ArgumentParser


def main(argv: list[str]) -> None:
    """Entry point for the CLI utility functions."""
    parser = ArgumentParser(
        prog="immunopipe utils",
        description="Utility commands to assist with immunopipe tasks"
    )
    #----------------------#
    # check-genes command  #
    #----------------------#
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
    #----------------------#
    # check-dim command    #
    #----------------------#
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
    #------------------------#
    # select-markers command #
    #------------------------#
    sm_command = parser.add_command(
        "select-markers",
        help=(
            "Select marker genes for each cluster from the ClusterMarkers process.\n"
            "This is useful for downstream analysis and visualization.\n"
        )
    )
    sm_command.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=str,
        help="Output directory for the pipeline.",
    )
    sm_command.add_argument(
        "--rscript",
        type=str,
        default="Rscript",
        help="Path to Rscript executable to run the R script",
    )
    sm_command.add_argument(
        "-t",
        "--top-n",
        type=int,
        default=10,
        help="Number of top marker genes to select for each cluster",
    )
    sm_command.add_argument(
        "--order-by",
        type=str,
        default="desc(avg_log2FC)",
        help=(
            "An expression to order the marker genes by.\n"
            "Available variables: avg_log2FC, pct.1, pct.2, p_val, p_val_adj.\n"
        ),
    )
    sm_command.add_argument(
        "-f",
        "--filter",
        type=str,
        default="p_val_adj < 0.05",
        help=(
            "An expression to filter the marker genes.\n"
            "Available variables: avg_log2FC, pct.1, pct.2, p_val, p_val_adj.\n"
            "Default: 'p_val_adj < 0.05'\n"
        ),
    )

    args = parser.parse_args(argv)
    if args.COMMAND == "check-genes":
        from .check_genes import check_genes
        check_genes(args)
    elif args.COMMAND == "check-dim":
        from .check_dim import check_dim
        check_dim(args)
    elif args.COMMAND == "select-markers":
        from .select_markers import select_markers
        select_markers(args)
