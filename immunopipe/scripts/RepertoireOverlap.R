library(immunarch)

immdata = "{{ in.immdata }}"
outdir = "{{ out.outdir }}"

dir.create(outdir, showWarnings = FALSE)
load(immdata)

imm_ov1 = repOverlap(immdata$data, .method="public", .verbose = F)
imm_ov2 = repOverlap(immdata$data, .method="jaccard", .verbose = F)

overlap_fig = file.path(outdir, 'repOverlap.overlap.png')
png(overlap_fig, width=1200, height=1000, res=100)
p = vis(imm_ov1)
print(p)
dev.off()

jaccard_fig = file.path(outdir, 'repOverlap.jaccard.png')
png(jaccard_fig, width=1200, height=1000, res=100)
p = vis(imm_ov2)
print(p)
dev.off()
