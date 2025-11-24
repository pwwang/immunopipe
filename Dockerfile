FROM justold/immunopipe-rpkgs:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER . /immunopipe

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /immunopipe
RUN fc-cache -f -v && \
    python -m pip install --no-cache-dir -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install --no-cache -v -E runinfo -E diagram -E dry && \
    python -m pip cache purge && \
    python -m poetry cache clear --all pypi && \
    pipen report update && \
    echo "cache=/tmp/npm-cache" > /home/mambauser/.npmrc && \
    ln -s /opt/conda/py_envs/numpy_v1/bin/python /opt/conda/bin/python_np1 && \
    ln -s $(micromamba run -n base python -c "import biopipen; print(biopipen.__path__[0])") \
        $(micromamba run -n base python_np1 -c "import sysconfig; print(sysconfig.get_path('purelib'))")

WORKDIR /workdir
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
