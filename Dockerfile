FROM justold/immunopipe-rpkgs:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /immunopipe

# scanpy needs numba, which for now requires numpy < 2.4
# see: https://github.com/numba/numba/issues/10263
RUN fc-cache -f -v && \
    python -m pip install --no-cache-dir -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install --no-cache -v --all-extras && \
    python -m pip install "numpy==2.3" && \
    python -m pip cache purge && \
    python -m poetry cache clear --all pypi && \
    rm -rf /home/$MAMBA_USER/.cache/pypoetry/* && \
    pipen report update && \
    rm -rf /home/$MAMBA_USER/.npm/* && \
    echo "cache=/tmp/npm-cache" > /home/mambauser/.npmrc

WORKDIR /workdir
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
