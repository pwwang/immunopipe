ARG BIOPIPEN_TAG=dev
FROM biopipen/base:${BIOPIPEN_TAG} AS builder

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /immunopipe

RUN micromamba env update -n base -f /immunopipe/docker/environment.yml && \
    python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v -E runinfo -E diagram && \
    pipen report update

# Final stage
FROM biopipen/base:${BIOPIPEN_TAG}

COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /immunopipe /immunopipe
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /opt/conda /opt/conda

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN fc-cache -f -v && \
    echo "cache=/tmp/npm-cache" > /home/mambauser/.npmrc


WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
