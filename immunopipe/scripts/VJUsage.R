library(dplyr)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)

samples = "{{ in.samples }}"
immdata = "{{ in.immdata }}"
outdir = "{{ out.outdir }}"
vdjtools_dir = file.path(outdir, 'vdj_tools')
vdjtools_patch = "{{ args.vdjtools_patch }}"
vdjtools = "{{ args.vdjtools }}"
ncores = {{ args.ncores }}

# ncores > 4 causing problems
ncores = min(ncores, 4)
registerDoParallel(ncores)

dir.create(outdir, showWarnings = FALSE)
dir.create(vdjtools_dir, showWarnings = FALSE)


load(samples)
load(immdata)

paths = samples %>%
    filter(Type == 'scTCR') %>%
    select(Sample, Path)
rawfiles = as.list(paths$Path)
names(rawfiles) = paths$Sample

foreach(sample=names(rawfiles)) %dopar% {
    print(paste("Handling", sample))
    rawfile = rawfiles[[sample]]
    rawdata = read.table(rawfile, sep=",", header=T, row.names=NULL) %>%
        as_tibble() %>%
        select(contig_id, reads)
    rawdata = as.data.frame(rawdata)
    rownames(rawdata) = sub("_contig_", "", rawdata$contig_id)
    rawdata = rawdata[, -1, drop=F]

    df = immdata$data[[sample]] %>%
        transmute(
            contig_id=ContigID,
            count=1, frequency=Proportion,
            CDR3nt=CDR3.nt,
            CDR3aa=CDR3.aa,
            V=V.name,
            D=D.name,
            J=J.name
        ) %>%
        tibble::as.tibble() %>%
        rowwise() %>%
        mutate(
            contig_id=paste(
                unlist(strsplit(contig_id, ";"))[1:(str_count(CDR3nt, ";")+1)],
                collapse=";"
            )
        ) %>%
        separate_rows(contig_id, CDR3nt, CDR3aa, V, D, J, sep=";") %>%
        mutate(count=rawdata[contig_id, 'reads']) %>%
        dplyr::select(-1)

    vdjfile = paste0(vdjtools_dir, "/", sample, ".txt")
    write.table(
        df, vdjfile,
        sep="\t", quote=F, row.names = F, col.names=T
    )

    patient = samples %>%
        filter(Type == 'scTCR' & Sample == sample) %>%
        pull('Patient')

    patient_dir = file.path(outdir, patient)
    dir.create(patient_dir, showWarnings = FALSE)

    command = sprintf(
        "bash %s %s PlotFancyVJUsage --plot-type png %s %s/%s",
        vdjtools_patch, vdjtools, vdjfile, patient_dir, sample
    )
    system(command)
}
