import os
import pyreadr
import cmdy

from datar.all import f, group_split, group_by, nrow, distinct, pull, filter

samples = "{{ in.samples }}"
outdir = "{{ out.outdir }}"
perl = "{{ args.perl }}"
rscript = "{{ args.rscript }}"
master_loader = "{{ args.master_loader }}"
count_loader = "{{ args.count_loader }}"

os.makedirs(outdir, exist_ok=True)

perl = cmdy.perl.bake(_exe=perl)
rscript = cmdy.rscript.bake(_exe=rscript)

rdata = pyreadr.read_r(samples)

# pull the patients with both scTCR and scRNA-seq data
# and with both Sources (typically Cancer and Normal)
for i, sdata in enumerate(
        rdata['samples'] >> group_by(f.Patient) >> group_split()
):
    # Skip if it is not paired
    # Sample	Type	Patient	Source	Path
    # MM003BM-earlier	scRNA	MM003-earlier	BM	...
    # MM003BM-earlier	scTCR	MM003-earlier	BM	...
    # MM003WBC-earlier	scRNA	MM003-earlier	WBC	...
    # MM003WBC-earlier	scTCR	MM003-earlier	WBC	...
    if nrow(sdata) < 4:
        continue
    patient = sdata.Patient.values[0]
    print(f'Doing patient: {patient} ...')

    prefixes = sdata >> distinct(f.Prefix) >> pull(f.Prefix, to='list')
    vdj_samples = sdata >> filter(f.Type == 'scTCR') >> pull(f.Path, to='list')
    rna_samples = sdata >> filter(f.Type == 'scRNA') >> pull(f.Path, to='list')
    perl(
        master_loader,
        _=[*vdj_samples, *rna_samples],
        o=outdir,
        s=patient,
        p=','.join(prefixes)
    ).fg

    rscript(
        count_loader,
        patient,
        *[prefix.split('-')[0] for prefix in prefixes],
        outdir
    )
