FROM justold/immunopipe-rpkgs:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /immunopipe

RUN fc-cache -f -v && \
    python -m pip install --no-cache-dir -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install --no-cache -v --all-extras && \
    python -m pip cache purge && \
    python -m poetry cache clear --all pypi && \
    rm -rf /home/$MAMBA_USER/.cache/pypoetry/* && \
    pipen report update && \
    rm -rf /home/$MAMBA_USER/.npm/* && \
    echo "cache=/tmp/npm-cache" > /home/mambauser/.npmrc && \
    R -e "install.packages('stringfish', repos='https://cloud.r-project.org')"
    # stringfish required by r-qs2, errored by conda-forge/r's stringfish v0.17.0
    # (v0.18 can't be installed with r4.4)
    # stringfish.so: undefined symbol: _ZTIN3tbb4taskE

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
