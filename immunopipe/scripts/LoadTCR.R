library(dplyr)
library(immunarch)

samples = "{{ in.samples }}"
outfile = "{{ out.outfile }}"
outdir = normalizePath("{{ job.outdir }}") # to absolute

load(samples)
all_sample_dir = file.path(outdir, "all_samples")
unlink(all_sample_dir, recursive = TRUE)
dir.create(all_sample_dir)

paths = samples %>% filter(Type == 'scTCR') %>% pull('Path')
for (sample in paths) {
    destfile = file.path(all_sample_dir, basename(sample))
    file.symlink(sample, destfile)
}

immdata = repLoad(all_sample_dir)
names(immdata$meta) = "File"
immdata$meta = immdata$meta %>% left_join(
    select(samples, c(1, 4:ncol(samples))) %>% mutate(
        File=tools::file_path_sans_ext(basename(Path))
    )
) %>% select(!c("File", "Path"))
names(immdata$data) = immdata$meta$Sample

save(immdata, file=outfile)
