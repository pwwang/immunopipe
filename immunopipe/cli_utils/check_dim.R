check_dim <- function(outdir) {
    cell_qc_file <- file.path(outdir, "qc", "cell_qc.txt")
    gene_qc_file <- file.path(outdir, "qc", "gene_qc.txt")
    if (!file.exists(cell_qc_file) || !file.exists(gene_qc_file)) {
        stop("QC files not found in the specified output directory.")
    }
    cell_qc <- read.table(cell_qc_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    cell_qc <- tibble::as_tibble(cell_qc)
    gene_qc <- read.table(gene_qc_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    gene_qc <- tibble::as_tibble(gene_qc)
    cat("Cell QC dimensions:\n")
    # print full content
    print(cell_qc, n=Inf)
    cat("\nGene QC dimensions:\n")
    print(gene_qc, n=Inf)
}