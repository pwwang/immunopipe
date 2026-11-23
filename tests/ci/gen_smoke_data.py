#!/usr/bin/env python3
"""Generate a tiny synthetic 10x (CellRanger v3) dataset for the CI smoke test.

Pure standard library (no numpy) so it can run inside any freshly created
conda environment. Deterministic: same seed -> same data.

Layout produced under <outdir>:
    data/prepared/<sample>/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
    data/samples.txt                (Sample, RNAData, Condition)
    smoke.config.toml               (documented minimal configuration)
"""

from __future__ import annotations

import gzip
import random
import sys
from pathlib import Path

SEED = 20260923
N_SAMPLES = 4
CELLS_PER_SAMPLE = 50
GENES_BG = 300
MARKERS_PER_CLUSTER = 25
N_CLUSTERS = 3
MT_GENES = ["MT-ND1", "MT-CO1", "MT-ATP6", "MT-CYB"]
LAMBDA_BG = 0.35
LAMBDA_MARKER = 3.0
LAMBDA_MT = 1.0

SAMPLES = [
    ("C1", "Control"),
    ("C2", "Control"),
    ("T1", "Treated"),
    ("T2", "Treated"),
]


def poisson(rng: random.Random, lam: float) -> int:
    """Knuth's algorithm; lambdas here are small (<10) so this is fine."""
    limit = 2.718281828459045 ** (-lam)
    k, p = 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


def write_gz(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="\n") as fh:
        fh.write(text)


def main() -> int:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "smoke").resolve()
    rng = random.Random(SEED)

    genes = [f"GENE{i:04d}" for i in range(1, GENES_BG + 1)]
    markers = [
        [f"CLU{k}_MARK{i:02d}" for i in range(1, MARKERS_PER_CLUSTER + 1)]
        for k in range(N_CLUSTERS)
    ]
    features = genes + [g for grp in markers for g in grp] + MT_GENES
    n_genes = len(features)

    total_cells = 0
    for sample, condition in SAMPLES:
        sdir = outdir / "data" / "prepared" / sample
        barcodes = [f"{sample}_{i:04d}-1" for i in range(1, CELLS_PER_SAMPLE + 1)]
        entries: list[str] = []
        for cell_idx, _ in enumerate(barcodes):
            # every cell belongs to one of N_CLUSTERS expression programs
            cluster = cell_idx % N_CLUSTERS
            high = set(markers[cluster])
            for gene_idx, gene in enumerate(features, start=1):
                if gene in high:
                    lam = LAMBDA_MARKER
                elif gene in MT_GENES:
                    lam = LAMBDA_MT
                else:
                    lam = LAMBDA_BG
                count = poisson(rng, lam)
                if count:
                    # 10x mtx is gene(row) x cell(column), 1-based
                    entries.append(f"{gene_idx} {cell_idx + 1} {count}")
        total_cells += len(barcodes)
        write_gz(
            sdir / "matrix.mtx.gz",
            "%%MatrixMarket matrix coordinate integer general\n"
            f"%\n{n_genes} {len(barcodes)} {len(entries)}\n"
            + "\n".join(entries)
            + "\n",
        )
        write_gz(sdir / "barcodes.tsv.gz", "\n".join(barcodes) + "\n")
        write_gz(
            sdir / "features.tsv.gz",
            "".join(f"{g}\t{g}\tGene Expression\n" for g in features),
        )

    samples_txt = "Sample\tRNAData\tCondition\n" + "".join(
        f"{s}\tdata/prepared/{s}\t{cond}\n" for s, cond in SAMPLES
    )
    (outdir / "data").mkdir(parents=True, exist_ok=True)
    (outdir / "data" / "samples.txt").write_text(samples_txt)

    # A tiny local GMT keeps the marker-enrichment step (ClusterMarkers always
    # runs it) offline and deterministic - no Enrichr web service involved.
    gmt = outdir / "data" / "genesets.gmt"
    gmt.write_text(
        "".join(
            f"SYNTH_CLU{k}\tsynthetic marker set {k}\t"
            + "\t".join(grp)
            + "\n"
            for k, grp in enumerate(markers)
        )
    )

    (outdir / "smoke.config.toml").write_text(
        "\n".join(
            [
                "# Minimal configuration, same shape as the documented",
                "# ImmunopipeMinimal.config.toml, pointed at the synthetic data.",
                'name = "ImmunopipeSmoke"',
                'outdir = "smoke-output"',
                "",
                "[SampleInfo.in]",
                'infile = ["data/samples.txt"]',
                "",
                "# Keep the smoke run small: 2 cores, and a local GMT instead of the",
                "# default Enrichr databases so the run needs no web service.",
                "[SeuratPreparing.envs]",
                "ncores = 2",
                "",
                "[SeuratClustering.envs]",
                "ncores = 2",
                "",
                "[ClusterMarkers.envs]",
                "ncores = 2",
                f'dbs = ["{gmt}"]',
                "",
                "[SeuratClusterStats.envs]",
                "ncores = 2",
                "",
            ]
        )
    )

    print(
        f"generated {total_cells} cells x {n_genes} genes across "
        f"{N_SAMPLES} samples in {outdir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
