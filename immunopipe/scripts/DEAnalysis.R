library(Seurat)
library(dplyr)
library(RcppTOML)
library(future)
library(enrichR)

samples = "{{ in.samples }}"
exprdir = "{{ in.exprdir }}"
outdir = "{{ out.outdir }}"
config = '{{ args.config }}'
ncores = {{ args.ncores }}

setEnrichrSite("Enrichr")
plan(strategy = "multicore", workers = ncores)
dir.create(outdir, showWarnings = FALSE)
load(samples) # samples
config = parseTOML(config, fromFile=FALSE)

# load all matrix and create seurat object for each sample
data = samples %>%
    filter(Type == 'scRNA') %>%
    mutate(Mat = file.path(exprdir, paste0(Prefix, ".mat.rds"))) #%>%
    # filter(file.exists(Mat)) %>%
    # right_join(metadata, by = c("Sample", "Patient"), suffix = c(".x", ""))

# Load all counts to seurat object for each sample
sample_names = as.character(data$Sample)
seurats = lapply(data$Mat, function(mat) {
    counts = readRDS(mat)
    CreateSeuratObject(counts=counts, min.cells=3, min.features=200)
})
names(seurats) = sample_names

do_one_case = function(name) {
    casedir = file.path(outdir, name)
    dir.create(casedir, showWarnings = FALSE)

    # filter data with filters
    filters = config[[name]]$filter
    casename = names(config[[name]]$case)[1]
    caseval = config[[name]]$case[[casename]]

    one_data = data
    for (fname in names(filters)) {
        fval = filters[[fname]]
        one_data = one_data[one_data[[fname]] == fval,,drop=F]
    }

    # check if case still has >= 2 levels
    if (
        length(unique(one_data[[casename]])) < 2 ||
        !caseval %in% one_data[[casename]]
    ) {
        stop(paste("Case value doesn't exist or not enough levels to perform DE analysis."))
    }

    seurat_obj = NULL
    seurat_objs = list()
    for (i in seq_len(nrow(one_data))) {
        sobj = seurats[[as.character(unlist(one_data[i, 'Sample']))]]
        sobj$group = as.character(unlist(one_data[i, casename]))
        if (i==1) {
            seurat_obj = sobj
        } else {
            seurat_objs = c(seurat_objs, list(sobj))
        }
    }

    de = merge(seurat_obj, seurat_objs)
    de$percent.mt = PercentageFeatureSet(de, pattern = '^MT-')
    de$percent.mt = de$percent.mt[de$percent.mt < 7.5]
    de <- NormalizeData(de, normalization.method = "LogNormalize")
    Idents(de) = "group"
    markers = FindMarkers(object = de, ident.1 = caseval)

    sig_markers = markers %>% filter(p_val_adj < 0.05)
    write.table(
        sig_markers,
        file.path(casedir, 'markers.txt'),
        row.names=T,
        col.names=T,
        sep="\t",
        quote=F
    )

    genes = rownames(sig_markers)
    enriched = enrichr(genes, config[[name]]$dbs)

    for (db in config[[name]]$dbs) {
        outtable = file.path(casedir, paste0('enrichr_', db, '.txt'))
        outfig = file.path(casedir, paste0('enrichr_', db, '.png'))

        write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

        png(outfig, width=1000, height=1000, res=100)
        print(plotEnrich(enriched[[db]], title=db))
        dev.off()
    }

}

for (name in names(config)) {
    do_one_case(name)
}
