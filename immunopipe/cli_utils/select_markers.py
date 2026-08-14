import sys
from pathlib import Path
from argx import Namespace

from .utils import run_r


def select_markers(args: Namespace) -> None:
    """Select marker genes for each cluster from the ClusterMarkers process."""

    outdir = Path(args.outdir)
    markers_dirs = list(outdir.joinpath("ClusterMarkers").glob("*.markers/"))
    if len(markers_dirs) == 0:
        print(f"No markers directories found in {outdir}/ClusterMarkers")
        sys.exit(1)

    if len(markers_dirs) > 1:
        print(f"Multiple markers directories found in {outdir}/ClusterMarkers:")
        for d in markers_dirs:
            print(f"  {d}")
        print("Please specify the correct output directory with --outdir")
        sys.exit(1)

    case_dirs = list(markers_dirs[0].glob("*/"))
    if len(case_dirs) == 0:
        print(f"No case directories found in {markers_dirs[0]}")
        sys.exit(1)

    if len(case_dirs) > 1:
        print(f"Multiple case directories found in {markers_dirs[0]}:")
        for d in case_dirs:
            print(f"  {d}")
        print("Please specify the correct output directory with --outdir")
        sys.exit(1)

    marker_files = case_dirs[0].rglob("*/markers.tsv")
    script = """
        library(rlang)
        library(dplyr)

        select_markers <- function(marker_file, top_n=10, order_by="desc(p_val_adj)") {
            markers <- read.table(marker_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
            cluster_name <- colnames(markers)[ncol(markers)]
            markers %%>%%
                group_by(!!sym(cluster_name)) %%>%%
                arrange(!!parse_expr(order_by)) %%>%%
                slice_head(n=top_n) %%>%%
                ungroup()
        }

        marker_files <- c(%(marker_files)s)
        top_n <- %(top_n)d
        order_by <- "%(order_by)s"
        markers <- do.call(rbind, lapply(marker_files, select_markers, top_n=top_n, order_by=order_by))
        write.table(markers, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    """ % {
        "marker_files": ", ".join(f'"{f}"' for f in marker_files),
        "top_n": args.top_n,
        "order_by": args.order_by,
    }  # noqa: E501
    rc, stdout, stderr = run_r(args.rscript, script)
    if rc != 0:
        print("Error running R script:")
        if stderr:
            print(stderr)
        sys.exit(rc)
    if stdout:
        print(stdout)
