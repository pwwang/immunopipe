FROM justold/immunopipe-rpkgs:latest AS builder

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /immunopipe

RUN fc-cache -f -v && \
    python -m pip install --no-cache-dir -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install --no-cache -v --all-extras && \
    python -m pip cache purge && \
    python -m poetry cache clear --all pypi && \
    pipen report update && \
    # --- Post-Install Cleanup ---
    # Remove JS source maps
    find /opt/conda -follow -type f -name '*.js.map' -delete && \
    # Remove Cython source files (not needed at runtime)
    find /opt/conda/lib/python*/site-packages -follow -type f -name '*.pyx' -delete && \
    # Remove Python type stub files (not needed at runtime)
    find /opt/conda/lib/python*/site-packages -follow -type f -name '*.pyi' -delete && \
    # Remove __pycache__ directories (includes all .pyc files)
    find /opt/conda -name '__pycache__' -type d -exec rm -rf '{}' + 2>/dev/null || true && \
    # Remove Python test directories in site-packages
    find /opt/conda/lib/python*/site-packages -maxdepth 2 \
        \( -name 'tests' -o -name 'test' \) -type d -exec rm -rf '{}' + 2>/dev/null || true && \
    # --- npm / node_modules Cleanup (pipen-report frontend) ---
    # Remove markdown docs and TypeScript declarations (not needed at runtime)
    find /opt/conda/lib/python*/site-packages/pipen_report/frontend/node_modules \
        \( -name '*.md' -o -name '*.d.ts' \) -type f -delete 2>/dev/null || true && \
    # Remove TypeScript source files (pre-built bundles are kept, sources are not needed)
    find /opt/conda/lib/python*/site-packages/pipen_report/frontend/node_modules \
        -name '*.ts' -not -name '*.d.ts' -type f -delete 2>/dev/null || true && \
    # Remove Flow type annotation and CoffeeScript source files
    find /opt/conda/lib/python*/site-packages/pipen_report/frontend/node_modules \
        \( -name '*.flow' -o -name '*.coffee' \) -type f -delete 2>/dev/null || true && \
    # Remove test directories inside node_modules
    find /opt/conda/lib/python*/site-packages/pipen_report/frontend/node_modules \
        -type d \( -name 'test' -o -name 'tests' -o -name '__tests__' \) \
        -exec rm -rf '{}' + 2>/dev/null || true && \
    # Remove .github directories inside node_modules
    find /opt/conda/lib/python*/site-packages/pipen_report/frontend/node_modules \
        -type d -name '.github' -exec rm -rf '{}' + 2>/dev/null || true && \
    # Remove the sample input from celltypist
    rm -rf /opt/conda/lib/python*/site-packages/celltypist/data/samples/Ensembl11* && \
    rm -rf /opt/conda/lib/python*/site-packages/celltypist/data/samples/GENCODE* && \
    rm -rf /opt/conda/lib/python*/site-packages/celltypist/data/samples/sample_cell_by_gene.csv && \
    # Remove datar's data
    rm -rf /opt/conda/lib/python*/site-packages/datar/data/*.gz && \
    # We don't need tensorboard
    rm -rf /opt/conda/lib/python*/site-packages/tensorboard* && \
    rm -rf /opt/conda/lib/python*/site-packages/tensorboard_data_server* && \
    # Remove plotnine's data
    rm -rf /opt/conda/lib/python*/site-packages/plotnine/data/* && \
    # Remove test files
    rm -rf /opt/conda/etc/conda/test-files/*

# Fresh stage: carry over only the installed environment and runtime files,
# not the full repo (docs/, tests/, skills/, notebooks, etc.)
FROM justold/immunopipe-rpkgs:latest

# Installed Python packages (poetry install target), including R scripts and
# report templates shipped as package data under site-packages/immunopipe/
COPY --from=builder /opt/conda /opt/conda
# The Python package itself: poetry installs it as an editable reference to
# /immunopipe, so this directory must exist in the final image.
# Only the package module is copied — docs/, tests/, skills/, etc. are excluded.
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /immunopipe/immunopipe /immunopipe/immunopipe
# Runtime entry script
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /immunopipe/docker /immunopipe/docker

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Redirect npm cache to a writable tmp location
RUN echo "cache=/tmp/npm-cache" > /home/$MAMBA_USER/.npmrc

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
