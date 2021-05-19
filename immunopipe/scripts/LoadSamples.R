library(dplyr)

samplefile = "{{ in.samplefile }}"
datadir = "{{ args.datadir }}"
outfile = "{{ out.outfile }}"
outdir = "{{ job.outdir }}"

# Sample	Type	Patient ClinicGroup1 CLinicGrouop2 ...
samples = read.table(
    samplefile,
    header=TRUE,
    row.names=NULL,
    sep="\t",
    check.names=FALSE
)

# prepend datadir to samples if possible
if (!datadir %in% c("None", "")) {
    datadir = normalizePath(datadir)
    samples = samples %>% mutate(
        Path=as.character(Path),
        Path=if_else(startsWith(Path, "/"), Path, file.path(datadir, Path))
    )
}
save(samples, file=outfile)

# export sample table
write.table(
    samples %>% mutate(File=basename(Path)) %>% select(-ncol(samples)),
    file.path(outdir, 'samples.txt'),
    row.names=FALSE,
    col.names=TRUE,
    sep="\t",
    quote=FALSE
)
