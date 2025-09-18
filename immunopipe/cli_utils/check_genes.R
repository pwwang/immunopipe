library(Seurat)
library(biopipen.utils)

check_genes <- function(seurat_path, genes, assay=NULL) {
    obj <- read_obj(seurat_path)
    assay_data <- GetAssayData(obj, assay=assay)
    all_genes <- rownames(obj)
    missing_genes <- setdiff(genes, all_genes)
    cat("Checked genes:\n")
    cat(paste0("  ", paste(genes, collapse=", "), "\n"))
    if (length(missing_genes) == 0) {
        cat("  All genes are present in the data.\n")
    } else {
        cat("Missing genes:\n")
        cat(paste0("  ", paste(missing_genes, collapse=", "), "\n"))
    }
}
