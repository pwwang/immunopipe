
library(RcppTOML)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(enrichR)
library(foreach)
library(doParallel)

pccdir = "{{ in.pccdir }}"
exprdir = "{{ in.exprdir }}"
samples = "{{ in.samples }}"
outdir = "{{ out.outdir }}"
tclusters = '{{ args.tclusters }}'
commoncfg = '{{ args.commoncfg }}'
config = '{{ args.config }}'
ncores = {{ args.ncores }}

registerDoParallel(ncores)
tclusters = parseTOML(tclusters, fromFile=FALSE)
config = parseTOML(config, fromFile=FALSE)

idents.tcell.pal = unlist(tclusters$colors)
names(idents.tcell.pal) = unlist(tclusters$names)

dir.create(outdir, showWarnings = FALSE)

load(samples)
tcrdir = "{{ in.tcrdir }}"
load(file.path(pccdir, "global.clonotype.ident.RData"))
# load(file.path(pccdir, "global.counts.RData"))
setEnrichrSite("Enrichr")

gsea_dbs = parseTOML(commoncfg, fromFile=FALSE)$GSEA_DBs
# patients = setdiff(names(multipt_samples), "COMPARING")
# normal_ident = clonotype.normal.ident[!is.na(clonotype.normal.ident)]
tumor_ident = clonotype.tumor.ident[!is.na(clonotype.tumor.ident)]
clusters = levels(tumor_ident)

mysamples = samples %>% filter(Type=="scRNA")
all_exprs = foreach (
    i = seq_len(nrow(mysamples)),
    .final = function(x) setNames(x, mysamples$Prefix)
) %dopar% {
    row = mysamples[i, ]
    pheno = unlist(strsplit(row$Prefix, "-", fixed=T))[1]
    phenofile = file.path(tcrdir, paste0(row$Patient, ".pheno.", pheno))
    if (!file.exists(phenofile)) {
        return (NULL)
    }
    pheno1 = read.table(phenofile, header=TRUE, row.names=NULL, sep="\t", check.names=FALSE) %>%
        select(barcode, clonotype) %>%
        filter(!is.na(clonotype))

    if (!file.exists(file.path(exprdir, paste0(row$Prefix, ".mat.rds")))) {
        return (NULL)
    }
    readRDS(file.path(exprdir, paste0(row$Prefix, ".mat.rds"))) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("barcode") %>%
        right_join(pheno1) %>%
        select(-"barcode") %>%
        group_by(clonotype) %>%
        summarise(across(everything(), sum)) %>%
        column_to_rownames("clonotype") %>%
        t()
}



is.normal = function(src, controls = c("wbc", "blood", "normal", "control")) {
    tolower(src) %in% controls
}

find_markers = function(prefs1, prefs2, cluster, level1, level2, casedir) {
    cldir = file.path(casedir, cluster)

    exprs1 = do.call(cbind, unname(all_exprs[prefs1]))
    exprs2 = do.call(cbind, unname(all_exprs[prefs2]))

    cl_clones = names(tumor_ident[tumor_ident == cluster])
    clones1 = intersect(cl_clones, colnames(exprs1))
    clones2 = intersect(cl_clones, colnames(exprs2))
    if (length(clones1) == 0 || length(clones2) == 0) {
        return (NULL)
    }
    exprs1 = exprs1[, clones1]
    exprs2 = exprs2[, clones2]


    de = tryCatch({
        s1 = CreateSeuratObject(counts=exprs1, min.cells=3, min.features=200)
        s1$group = level1
        s2 = CreateSeuratObject(counts=exprs2, min.cells=3, min.features=200)
        s2$group = level2
        de = merge(s1, s2)
        de$percent.mt = PercentageFeatureSet(de, pattern = '^MT-')
        de$percent.mt = de$percent.mt[de$percent.mt < 7.5]
        de <- NormalizeData(de, normalization.method = "LogNormalize")
        Idents(de) = "group"
        de
    }, error=function(e) {
        NULL
    })
    if (is.null(de)) { return(NULL) }

    markers = FindMarkers(object = de, ident.1 = level1)
    sig_markers = markers %>% filter(p_val_adj < 0.05)
    if (nrow(sig_markers) > 0) {
        dir.create(cldir, showWarnings=FALSE)

        write.table(
            sig_markers,
            file.path(cldir, 'markers.txt'),
            row.names=T,
            col.names=T,
            sep="\t",
            quote=F
        )

        genes = rownames(sig_markers)
        enriched = enrichr(genes, gsea_dbs)
        for (db in gsea_dbs) {
            outtable = file.path(cldir, paste0('enrichr_', db, '.txt'))
            outfig = file.path(cldir, paste0('enrichr_', db, '.png'))

            write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

            png(outfig, width=1000, height=1000, res=100)
            print(plotEnrich(enriched[[db]], title=db))
            dev.off()
        }
    }
}


do_one_case = function(name) {
    cfg = config[[name]]
    casedir = file.path(outdir, name)

    mysamples = samples %>% filter(Type=="scRNA")
    for(filtname in names(cfg$filter)) {
        mysamples = mysamples %>% filter(mysamples[[filtname]] == cfg$filter[[filtname]])
    }

    casename = names(cfg$case)[1]
    clinics = as.character(mysamples[[casename]])
    clinics = unique(clinics[!endsWith(clinics, "-U")])
    if (length(clinics) != 2) {
        warning(paste0('Warning: ', casename, ' does not have 2 levels exactly, skipping'))
        return (NULL)
    }

    dir.create(casedir, showWarnings = FALSE)
    level1 = clinics[1]
    level2 = clinics[2]

    prefs1 = mysamples %>% filter(mysamples[[casename]] == level1) %>% pull(Prefix)
    prefs2 = mysamples %>% filter(mysamples[[casename]] == level2) %>% pull(Prefix)

    for (cluster in clusters) {
        find_markers(prefs1, prefs2, cluster, level1, level2, casedir)
    }

}


for (name in names(config)) {
    do_one_case(name)
}
