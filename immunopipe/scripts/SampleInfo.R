library(dplyr)

samplefile = "{{ in.samplefile }}"
metafile = "{{ in.metafile }}"
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
if ('End' %in% names(samples)) {
    samples = samples %>% select(-End)
}

# Sample	Patient	...
metadata = read.table(
    metafile,
    header=TRUE,
    row.names=NULL,
    sep="\t",
    check.names=FALSE
)

metagroups = list()
metanames = colnames(metadata)[1:2]
for (idx in 3:ncol(metadata)) {
    name = colnames(metadata)[idx]
    if (grepl("\\\\.\\\\d+$", name)) {
        mname = gsub("\\\\.\\\\d+$", "", name)
        gidx = substring(name, nchar(mname)+2)
        if (gidx %in% names(metagroups)) {
            metagroups[[gidx]] = c(metagroups[[gidx]], mname)
        } else {
            metagroups[[gidx]] = mname
        }
        metanames = c(metanames, mname)
    } else {
        metagroups[[idx+1000]] = name
        metanames = c(metanames, name)
    }
}
colnames(metadata) = metanames

# prepend datadir to samples if possible
if (!datadir %in% c("None", "")) {
    datadir = normalizePath(datadir)
    samples = samples %>% mutate(
        Path=as.character(Path),
        Path=if_else(startsWith(Path, "/"), Path, file.path(datadir, Path)),
        Genes=as.character(Genes),
        Genes=if_else(
            startsWith(Genes, "/"),
            Genes,
            if (nchar(Genes) == 0) "" else file.path(datadir, Genes)
        ),
        Matrix=as.character(Matrix),
        Matrix=if_else(
            startsWith(Matrix, "/"),
            Matrix,
            if (nchar(Matrix) == 0) "" else file.path(datadir, Matrix)
        )
    )
}

# check if files exist
for (path in samples$Path) {
    if (!file.exists(path)) {
        stop(paste(path, "does not exist!"))
    }
}

# check if files exist
for (genes in samples$Genes) {
    if (nchar(genes) > 0 && !file.exists(genes)) {
        stop(paste(genes, "does not exist!"))
    }
}

# check if files exist
for (mat in samples$Matrix) {
    if (nchar(mat) > 0 && !file.exists(mat)) {
        stop(paste(mat, "does not exist!"))
    }
}

# move tcr files into one directory for immunarch to load
tcrdir = file.path(outdir, "tcr-files")
unlink(tcrdir, recursive = TRUE, force = TRUE)
dir.create(tcrdir, showWarnings = FALSE)
new_paths = c()
for (i in seq_len(nrow(samples))) {
    path = samples[i, 'Path']
    sample = as.character(samples[i, 'Sample'])
    parts = unlist(strsplit(basename(path), '.', fixed=TRUE))
    # 10x specific
    if (grepl("annotations", parts[1])) {
        parts = c(sample, parts)
    } else {
        parts[1] = sample
    }
    new_path = file.path(tcrdir, paste(parts, collapse='.'))
    file.symlink(path, new_path)
    new_paths = c(new_paths, new_path)
}
samples$Path = new_paths

# Add prefixes
samples = samples %>%
    group_by(Patient) %>%
    mutate(Prefix=paste(Source, cur_group_id(), sep='-'))

save(samples, metadata, tcrdir, metagroups, file=outfile)

# export sample table
write.table(
    samples %>% mutate(Path=basename(Path)) %>% select(1:5),
    file.path(outdir, 'samples.txt'),
    row.names=FALSE,
    col.names=TRUE,
    sep="\t",
    quote=FALSE
)

write.table(
    metadata,
    file.path(outdir, 'metadata.txt'),
    row.names=FALSE,
    col.names=TRUE,
    sep="\t",
    quote=FALSE
)
