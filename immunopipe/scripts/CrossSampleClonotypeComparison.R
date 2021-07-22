library(dplyr)
library(immunarch)
library(RcppTOML)
library(ggprism)
library(doParallel)

immdata = "{{ in.immdata }}"
outdir = "{{ out.outdir }}"
config = '{{ args.config }}'
ncores = {{ args.ncores }}

dir.create(outdir, showWarnings = FALSE)
registerDoParallel(ncores)
load(immdata) # samples, metadata

# [PatientSamples]
# MM005 = ["MM005-earlier", "MM005-postr"]
# MM006 = ["MM006-pre", "MM006-postr"]
#     [PatientSamples.COMPARING]
#     ORDER = ["Type", "Status"]
#     Type = ["BM", "WBC"]
#     Status = ["Pre", "Post"]
config = parseTOML(config, fromFile=FALSE)


clonotypes_composition = function(patient, patdir, pdata, meta, type) {
    pr = pubRep(pdata, type, .verbose = F)
    clonotypes = pr %>% rowwise() %>%
        mutate(rowsum=sum(c_across(-c(1:2)))) %>%
        arrange(desc(rowsum))

    write.table(
        clonotypes,
        file.path(patdir, paste('clonotypes', type, 'txt', sep='.')),
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )

    tc = trackClonotypes(
        pdata,
        clonotypes[[paste("CDR3", type, sep=".")]] %>% head(10),
        .col = type
    )
    p = vis(tc) +
        ggtitle(paste0(patient, " Top 10 shared clonotype composition (CDR3.", type, ")"))

    png(file.path(patdir, paste('clonotypes', type, 'png', sep='.')),
        width=1800, height=1000, res=100)
    print(p)
    dev.off()
}

gain_and_loss = function(data, sample1, sample2, suffices) {
    data1 = data[[sample1]][,c("CDR3.nt", "Proportion"),drop=F]
    data2 = data[[sample2]][,c("CDR3.nt", "Proportion"),drop=F]
    ret = full_join(data1, data2, by="CDR3.nt", suffix = paste0("_", suffices))
    prop1 = ret[[paste("Proportion", suffices[1], sep='_')]]
    prop2 = ret[[paste("Proportion", suffices[2], sep='_')]]

    ret = ret %>%
        mutate(
            FreqChange=if_else(
                is.na(prop2),
                -prop1,
                if_else(
                    is.na(prop1),
                    prop2,
                    prop2 - prop1
                )
            ),
            Type=if_else(
                is.na(prop2),
                "Vanish",
                if_else(
                    is.na(prop1),
                    "Emerge",
                    if_else(FreqChange>0, "Expand", "Collapse")
                )
            )
        ) %>%
        arrange(desc(abs(FreqChange)))
    cret = ret %>% group_by(Type) %>% count()

    list(data=ret, count=cret)
}


plot_gain_loss = function(gldata, top) {
    ggplot(gldata, aes(x=reorder(CDR3.nt, abs(FreqChange)), y=FreqChange, fill=Type)) +
        geom_col() +
        theme_prism(base_size=16) +
        scale_fill_prism("colors") +
        coord_flip() +
        ggtitle(paste("Top", top, "vanished or emerged clonotypes"))
}

for (patient in names(config))  {
    if (patient == "COMPARING") next
    patdir = file.path(outdir, patient)
    dir.create(patdir, showWarnings = FALSE)

    meta = immdata$meta %>%
        filter(Patient %in% config[[patient]])

    metanames = c()
    for (metaname in names(config$COMPARING)) {
        if (metaname %in% c('ORDER', 'DE_DBS', 'DE_TOP_CLONOTYPES')) next
        meta[[metaname]] = factor(
            meta[[metaname]],
            levels = config$COMPARING[[metaname]]
        )
        metanames = c(metanames, metaname)
    }
    meta = arrange(meta, meta[, metanames])

    # nt
    pdata = immdata$data[meta$Sample]
    clonotypes_composition(patient, patdir, pdata, meta, 'nt')
    clonotypes_composition(patient, patdir, pdata, meta, 'aa')

    # Clonotype changes
    compairing_metas = config$COMPARING$ORDER
    for (m1 in seq(1, length(compairing_metas)-1)) {
        mname1 = compairing_metas[m1]
        for (mval in config$COMPARING[[mname1]]) {
            meta_cc = meta %>% filter(meta[[mname1]] == mval)

            for (m2 in seq(m1+1, length(compairing_metas))) {
                mname2 = compairing_metas[m2]
                comps = config$COMPARING[[mname2]]
                msamples = meta_cc %>%
                    arrange(meta_cc[[mname2]]) %>%
                    pull(Sample)
                gldata = gain_and_loss(pdata, msamples[1], msamples[2], comps)

                glfile = file.path(
                    patdir,
                    paste0("gainloss_", mval, '_', paste(comps, collapse="vs"), ".txt")
                )
                glcfile = file.path(
                    patdir,
                    paste0("glcount_", mval, '_', paste(comps, collapse="vs"), ".txt")
                )
                write.table(gldata$data, glfile, col.names=T, row.names= F, sep="\t", quote=F)
                write.table(gldata$count, glcfile, col.names=T, row.names= F, sep="\t", quote=F)

                glfig = file.path(
                    patdir,
                    paste0("gainloss_", mval, '_', paste(comps, collapse="vs"), ".png")
                )

                g = plot_gain_loss(
                    gldata$data %>%
                    filter(Type %in% c("Vanish", "Emerge")) %>%
                    head(30) %>%
                    mutate(CDR3.nt = paste(substr(CDR3.nt, 1, 30), '...')),
                    top=30
                )
                png(glfig, width=1200, height=1000, res=100)
                print(g)
                dev.off()
            }
        }
    }
}
