library(dplyr)
library(immunarch)

samples = "{{ in.samples }}"
outfile = "{{ out.outfile }}"
outdir = normalizePath("{{ job.outdir }}") # to absolute

load(samples)

immdata = repLoad(tcrdir)
immdata$meta$Sample = gsub(
    '.filtered_contig_annotations',
    '',
    immdata$meta$Sample,
    fixed=TRUE
)
names(immdata$data) = immdata$meta$Sample
# add Source to immdata$meta
immdata$meta = left_join(
        immdata$meta,
        samples %>% filter(Type == 'scTCR')
    )

save(immdata, file=outfile)
