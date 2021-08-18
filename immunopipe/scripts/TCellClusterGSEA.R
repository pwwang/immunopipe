library(dplyr)
library(enrichR)

ctdir = "{{ in.ctdir }}" # T cell cluster directory
outdir = "{{ out.outdir }}"

setEnrichrSite("Enrichr")

dir.create(outdir, showWarnings = FALSE)

DBS = c("KEGG_2021_Human")

clusterGSEA = function(rdfile) {
    ident = unlist(strsplit(basename(rdfile), ".", fixed=T))[2]
    load(rdfile)

    genes = rownames(markers)

    for (db in DBS) {
        enriched = enrichr(genes, db)
        outtable = file.path(outdir, paste0('Ident.',ident,'.', db, '.txt'))
        outfig = file.path(outdir, paste0('Ident.',ident,'.', db, '.png'))

        write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

        png(outfig, width=1000, height=1000, res=100)
        print(plotEnrich(enriched[[db]], title=db))
        dev.off()
    }
}

for (rdfile in Sys.glob(file.path(ctdir, "Markers", "*.markers.RData"))) {
    clusterGSEA(rdfile)
}
