library(Seurat)
library(dplyr)
library(enrichR)
library(RcppTOML)
library(foreach)
library(doParallel)

samples = "{{ in.samples }}"
immdata = "{{ in.immdata }}"
exprdir = "{{ in.exprdir }}"
ccdir = "{{ in.ccdir }}"
outdir = "{{ out.outdir }}"
config = '{{ args.config }}'
ncores = {{ args.ncores }}

registerDoParallel(ncores)
setEnrichrSite("Enrichr")
dir.create(outdir, showWarnings = FALSE)
load(immdata) # samples, metadata
load(samples) # samples, metadata

# [PatientSamples]
# MM005 = ["MM005-earlier", "MM005-postr"]
# MM006 = ["MM006-pre", "MM006-postr"]
#     [PatientSamples.COMPARING]
#     ORDER = ["Type", "Status"]
#     Type = ["BM", "WBC"]
#     Status = ["Pre", "Post"]
config = parseTOML(config, fromFile=FALSE)
top = config$COMPARING$DE_TOP_CLONOTYPES

create_seurat_objs = function(prefix, clones) {
    exprs = readRDS(file.path(exprdir, paste0(prefix, ".mat.rds")))
    exprs1 = exprs[, intersect(colnames(exprs), paste(prefix, clones, sep="_"))]
    exprs2 = exprs[, setdiff(colnames(exprs), colnames(exprs1))]
    s1 = CreateSeuratObject(counts=exprs1, min.cells=3, min.features=200)
    s1$group = 'with'
    s2 = CreateSeuratObject(counts=exprs2, min.cells=3, min.features=200)
    s2$group = 'without'
    merge(s1, s2)
}

get_clonotypes = function(sample, clone) {
    immdata$data[[sample]] %>%
        filter(CDR3.nt == clone) %>%
        pull(Barcode) %>%
        strsplit(';', fixed=TRUE) %>%
        unlist()
}

get_prefix_clones = function(pat, type, comp, sample_prefs) {
    ccfile = file.path(
        ccdir,
        pat,
        paste0(paste('gainloss', type, paste(comp, collapse='vs'), sep='_'), '.txt')
    )

    clonedata = read.table(
        ccfile,
        header=TRUE,
        row.names=NULL,
        sep="\t",
        check.names=FALSE
    ) %>%
        filter(Type %in% c('Vanish', 'Emerge')) %>%
        slice_head(n=top)
    # sample_prefs
    # MM001BM-earlier => BM-1
    # MM001BM-post => BM-2

    sams = names(sample_prefs)

    out = list()
    for (i in seq_len(nrow(clonedata))) {
        clone = as.character(clonedata[i, 1])
        prop1 = clonedata[i, 2]
        prop2 = clonedata[i, 3]
        sam = if (is.na(prop1)) sams[2] else sams[1]
        pref = sample_prefs[[sam]]
        out = c(out, list(list(
            pref = pref,
            clone = clone,
            sams = sams,
            props = c(prop1, prop2),
            clonotypes = get_clonotypes(sam, clone)
        )))
    }
    out
}

get_samples_by_patient = function(patients, type) {
    filtered = samples %>%
        filter(Patient %in% patients, Source == type)
    out = filtered %>% pull(Prefix) %>% unique()
    outnames = filtered %>% pull(Sample) %>% unique() %>% unlist() %>% as.character()
    names(out) = outnames
    out
}

do_pref = function(pref, casedir, clone, clonotypes, sams, props) {
    clonedir = file.path(casedir, clone)
    dir.create(clonedir, showWarnings = FALSE)

    infofile = file.path(clonedir, 'info.txt')
    write.table(
        data.frame(Sample=sams, Proportion=props),
        infofile,
        row.names=FALSE,
        col.names=TRUE,
        sep="\t",
        quote=FALSE
    )

    de = create_seurat_objs(pref, clonotypes)
    de$percent.mt = PercentageFeatureSet(de, pattern = '^MT-')
    de$percent.mt = de$percent.mt[de$percent.mt < 7.5]
    de <- NormalizeData(de, normalization.method = "LogNormalize")
    Idents(de) = "group"

    markers = FindMarkers(object = de, ident.1 = 'with')

    sig_markers = markers %>% filter(p_val_adj < 0.05)
    write.table(
        sig_markers,
        file.path(clonedir, 'markers.txt'),
        row.names=T,
        col.names=T,
        sep="\t",
        quote=F
    )

    genes = rownames(sig_markers)
    enriched = enrichr(genes, config$COMPARING$DE_DBS)

    for (db in config$COMPARING$DE_DBS) {
        outtable = file.path(clonedir, paste0('enrichr_', db, '.txt'))
        outfig = file.path(clonedir, paste0('enrichr_', db, '.png'))

        write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

        png(outfig, width=1000, height=1000, res=100)
        print(plotEnrich(enriched[[db]], title=db))
        dev.off()
    }
}

do_one_patient_one_type = function(pat, patients, type, comp) {
    patdir = file.path(outdir, pat)
    dir.create(patdir, showWarnings = FALSE)
    casedir = file.path(patdir, paste(type, paste(comp, collapse='vs'), sep="_"))
    dir.create(casedir, showWarnings = FALSE)

    sample_prefs = get_samples_by_patient(patients, type)
    # MM001BM-earlier => BM-1
    # MM001BM-post => BM-2
    pref_clones = get_prefix_clones(pat, type, comp, sample_prefs)
    # pref => list(
    #    clone=TGTGCTGTGGTGGACCCCCCTAACCAGGCAGGAACTGCTCTGATCTTT;TGTGCCAGCGAAAAACTTTCCTACGAGCAGTACTTC
    #    clonotypes=c('TGTGCTGTGGTGGACCCC-1', 'GGAACTGCTCTGATCTTT-1', ...)
    # )

    foreach (i=seq_along(pref_clones)) %dopar% {
        prefc = pref_clones[[i]]
        print(paste('  * Do clone', prefc$pref, prefc$clone))
        do_pref(
            prefc$pref,
            casedir,
            prefc$clone,
            prefc$clonotypes,
            prefc$sams,
            prefc$props
        )
    }
}

for (pat in names(config))  {
    if (pat == "COMPARING") next
    compname = config$COMPARING$ORDER[2]

    print(paste('Got patient', pat, '...'))
    for (type in config$COMPARING$Source) {
        print(paste('- Do type', type, '...'))
        do_one_patient_one_type(
            pat, # MM005
            config[[pat]], # MM005-earlier, # MM005-postr
            type,
            config$COMPARING[[compname]]
        )
    }
}
