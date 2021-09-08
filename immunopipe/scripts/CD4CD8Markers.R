
library(RcppTOML)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(enrichR)
library(foreach)
library(doParallel)

cldir = "{{ in.cldir }}"
outdir = "{{ out.outdir }}"
CD4CD8clusters = '{{ args.cd4cd8clusters }}'
commoncfg = '{{ args.commoncfg }}'
ncores = {{ args.ncores }}

registerDoParallel(ncores)
CD4CD8clusters = parseTOML(CD4CD8clusters, fromFile=FALSE)
dir.create(outdir, showWarnings = FALSE)

# load(file.path(pccdir, "global.counts.RData"))
setEnrichrSite("Enrichr")

gsea_dbs = parseTOML(commoncfg, fromFile=FALSE)$GSEA_DBs
load(file.path(cldir, "Seurat-1", 'global.tcell.RData'))

find_markers = function(name, clusters, casedir) {
    markers <- FindMarkers(global.obj, ident.1=clusters)
    save(
        markers,
        file=file.path(casedir, paste0("markers.RData"))
    )
    write.table(
        markers,
        file.path(casedir, paste0("markers.txt")),
        col.names = TRUE,
        row.names = TRUE,
        sep="\t",
        quote = FALSE
    )

    # GSEA
    dir.create(casedir, showWarnings = FALSE)
    genes = markers %>% filter(p_val_adj < 0.05) %>% rownames()
    enriched = enrichr(genes, gsea_dbs)

    for (db in gsea_dbs) {
        outtable = file.path(casedir, paste0('enrichr_', db, '.txt'))
        outfig = file.path(casedir, paste0('enrichr_', db, '.png'))

        write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

        png(outfig, width=1000, height=1000, res=100)
        print(plotEnrich(enriched[[db]], title=db))
        dev.off()
    }
}


do_one_case = function(name) {
    casedir = file.path(outdir, name)
    dir.create(casedir, showWarnings = FALSE)

    clusters = CD4CD8clusters[[name]]
    find_markers(name, clusters, casedir)
}


for (name in names(CD4CD8clusters)) {
    do_one_case(name)
}
