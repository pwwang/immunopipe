library(biopipen.utils)


infile <- {{in.infile | quote}}
outfile <- {{out.outfile | quote}}
sample_col <- {{envs.sample | quote}}

obj <- read_obj(infile)
if ("Sample" %in% colnames(obj@meta.data) && sample_col == "Sample") {
    # do nothing, make a symlink to outfile
    if (file.exists(outfile)) {
        file.remove(outfile)
    }
    file.symlink(infile, outfile)
} else if (!sample_col %in% colnames(obj@meta.data)) {
    stop(paste0("Sample column '", sample_col, "' not found in metadata"))
} else {
    obj@meta.data$Sample <- obj@meta.data[[sample_col]]
    save_obj(obj, outfile)
}
