FROM mambaorg/micromamba:1.4.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

# Install dependencies and make /data for future mounting
RUN micromamba install -y -n base -f /immunopipe/docker/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /data
WORKDIR /immunopipe
RUN fc-cache -f -v && \
    python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v -E runinfo -E diagram && \
    pipen report update && \
    echo "cache=/tmp/npm-cache" > /home/mambauser/.npmrc

WORKDIR /workdir
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
