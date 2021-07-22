library(immunarch)
library(RcppTOML)

samples = "{{ in.samples }}"
immdata = "{{ in.immdata }}"
outdir = "{{ out.outdir }}"
config = '{{ args.config }}'

dir.create(outdir, showWarnings = FALSE)
grouping = parseTOML(config, fromFile=FALSE)$grouping

load(samples)
load(immdata)

# Number of clonotypes of each sample
exp_vol = repExplore(immdata$data, .method = "volume")
noc_dir = file.path(outdir, 'Number_of_clonotypes')
dir.create(noc_dir, showWarnings = FALSE)

## Overall
ov_noc_dir = file.path(noc_dir, 'Overall')
dir.create(ov_noc_dir, showWarnings = FALSE)

overallfig = file.path(ov_noc_dir, 'Number_of_clonotypes.png')
png(overallfig, width=1200, height=1000, res=100)
p = vis(exp_vol)
print(p)
dev.off()

# sourcefig = file.path(ov_noc_dir, 'Number_of_clonotypes-Source.png')
# png(overallfig, width=1200, height=1000, res=100)
# p = vis(exp_vol, .by='Source', .meta=immdata$meta)
# print(p)
# dev.off()

# todo: 1-group, 2-groups, etc...
for (groups in grouping) {
    groupfig = file.path(
        noc_dir,
        paste0('Number_of_clonotypes-',paste(groups, collapse='.vs.'),'.png')
    )
    png(groupfig, width=1200, height=1000, res=100)
    p = vis(exp_vol, .by=groups, .meta=immdata$meta)
    print(p)
    dev.off()
}

# Distribution of clonotype abundance
exp_cnt <- repExplore(immdata$data, .method = "count")
doc_dir = file.path(outdir, 'Distribution_of_clonotype_abundance')
dir.create(doc_dir, showWarnings = FALSE)

fig = file.path(doc_dir, 'Distribution_of_clonotype_abundance.png')
png(fig, width=1000, height=1000, res=100)
p = vis(exp_cnt)
print(p)
dev.off()

# Top clones
imm_top = repClonality(immdata$data, .method='top', .head=c(10, 100, 1000, 3000))
top_dir = file.path(outdir, 'Top_clonal_proportion')
dir.create(top_dir, showWarnings = FALSE)

fig = file.path(top_dir, 'Top_clonal_proportion.png')
png(fig, width = 1200, height = 1000, res = 100)
p = vis(imm_top, .meta=immdata$meta)
print(p)
dev.off()

for (groups in grouping) {
    groupfig = file.path(
        top_dir,
        paste0('Top_clonal_proportion-',paste(groups, collapse='.vs.'),'.png')
    )
    png(groupfig, width=1200, height=1000, res=100)
    p = vis(imm_top, .by=groups, .meta=immdata$meta)
    print(p)
    dev.off()
}

# Rare clones
imm_rare = repClonality(immdata$data, .method='rare')
rare_dir = file.path(outdir, 'Rare_clonal_proportion')
dir.create(rare_dir, showWarnings = FALSE)

fig = file.path(rare_dir, 'Rare_clonal_proportion.png')
png(fig, width = 1200, height = 1000, res = 100)
p = vis(imm_rare, .meta=immdata$meta)
print(p)
dev.off()

for (groups in grouping) {
    groupfig = file.path(
        rare_dir,
        paste0('Rare_clonal_proportion-',paste(groups, collapse='.vs.'),'.png')
    )
    png(groupfig, width=1200, height=1000, res=100)
    p = vis(imm_rare, .by=groups, .meta=immdata$meta)
    print(p)
    dev.off()
}


# clonal space homeostasis
imm_hom = repClonality(immdata$data, .method='homeo')
hom_dir = file.path(outdir, 'Relative_clonal_abundance')
dir.create(hom_dir, showWarnings = FALSE)

fig = file.path(hom_dir, 'Relative_clonal_abundance.png')
png(fig, width = 1200, height = 1000, res = 100)
p = vis(imm_hom, .meta=immdata$meta)
print(p)
dev.off()

for (groups in grouping) {
    groupfig = file.path(
        hom_dir,
        paste0('Relative_clonal_abundance-',paste(groups, collapse='.vs.'),'.png')
    )

    png(groupfig, width=1200, height=1000, res=100)
    p = vis(imm_hom, .by=groups, .meta=immdata$meta)
    print(p)
    dev.off()
}
