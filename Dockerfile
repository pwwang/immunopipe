FROM justold/immunopipe-rpkgs:latest AS builder

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /immunopipe

RUN fc-cache -f -v && \
    uv pip install --system -e .[diagram,runinfo,dry,email,cli-gbatch] && \
    pipen report update && \
    python /immunopipe/docker/cleanup.py

# Fresh stage: carry over only the installed environment and runtime files,
# not the full repo (docs/, tests/, skills/, notebooks, etc.)
FROM justold/immunopipe-rpkgs:latest

# Installed Python packages (uv pip install target), including R scripts and
# report templates shipped as package data under site-packages/immunopipe/
COPY --from=builder /opt/conda /opt/conda
# The Python package itself: uv installs it as an editable reference to
# /immunopipe, so this directory must exist in the final image.
# Only the package module is copied — docs/, tests/, skills/, etc. are excluded.
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /immunopipe/immunopipe /immunopipe/immunopipe
# Runtime entry script
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /immunopipe/docker /immunopipe/docker

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Backport scverse/scanpy#4480 so scanpy works with umap-learn >= 0.6: 0.6.0 stopped
# copying the knn arrays, which made pp.neighbors fail on >= 8192 cells with a numba
# TypingError (lmcinnes/umap#1335) and killed CellTypeAnnotation's celltypist run.
# This lives in the final stage because scanpy is installed by the builder stage
# above. The script is idempotent and fails loudly if scanpy's code moved --
# drop this step once the image ships scanpy >= 1.12.6.
RUN python /immunopipe/docker/patch_scanpy_umap06.py

# Redirect npm cache to a writable tmp location
RUN echo "cache=/tmp/npm-cache" > /home/$MAMBA_USER/.npmrc

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
