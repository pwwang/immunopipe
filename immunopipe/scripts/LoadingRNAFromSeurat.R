library(rlang)
library(dplyr)
library(tidyseurat)
library(biopipen.utils)

infile <- {{in.infile | quote}}
outfile <- {{out.outfile | quote}}
sample_col <- {{envs.sample | quote}}
mutaters <- {{envs.mutaters | r}}
subset <- {{envs.subset | r}}
ncores <- {{envs.ncores | r}}

qs2::qopt("nthreads", value = ncores)

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

if (
    "Sample" %in% colnames(obj@meta.data) &&
    identical(sample_col, "Sample") &&
    (is.null(mutaters) || length(mutaters) == 0) &&
    is.null(subset)
) {
    # do nothing, make a symlink to outfile
    if (file.exists(outfile)) {
        file.remove(outfile)
    }
    file.symlink(infile, outfile)
} else {
    if (length(sample_col) > 1) {
        obj@meta.data$Sample <- apply(obj@meta.data[, sample_col, drop = FALSE], 1, paste, collapse = "_")
    } else {
        obj@meta.data$Sample <- obj@meta.data[[sample_col]]
    }
    obj <- MutateSeuratMeta(obj, mutaters)
    if (!is.null(subset)) {
       obj <- filter(obj, !!parse_expr(subset))
    }
    save_obj(obj, outfile)
}
