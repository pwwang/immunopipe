library(biopipen.utils)


infile <- {{in.infile | quote}}
outfile <- {{out.outfile | quote}}
sample_col <- {{envs.sample | quote}}

obj <- read_obj(infile)
if (!all(sample_col %in% colnames(obj@meta.data))) {
    missing_cols <- setdiff(sample_col, colnames(obj@meta.data))
    stop(
        paste0(
            "Sample column(s) '", paste(missing_cols, collapse = ", "), "' not found in metadata. Available columns: ",
            paste(colnames(obj@meta.data), collapse = ", ")
        )
    )
}

if ("Sample" %in% colnames(obj@meta.data) && identical(sample_col, "Sample")) {
    # do nothing, make a symlink to outfile
    if (file.exists(outfile)) {
        file.remove(outfile)
    }
    file.symlink(infile, outfile)
} else if (length(sample_col) > 1) {
    obj@meta.data$Sample <- apply(obj@meta.data[, sample_col, drop = FALSE], 1, paste, collapse = "_")
    save_obj(obj, outfile)
} else {
    obj@meta.data$Sample <- obj@meta.data[[sample_col]]
    save_obj(obj, outfile)
}
