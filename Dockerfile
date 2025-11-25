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
    echo -e '#!/bin/bash\nexec /opt/conda/py_envs/numpy_v1/bin/python "$@"' > /opt/conda/bin/python_np1 && \
    chmod +x /opt/conda/bin/python_np1 && \
    ln -s $(/opt/conda/bin/python -c "import biopipen; print(biopipen.__path__[0])") \
        $(/opt/conda/bin/python_np1 -c "import sysconfig; print(sysconfig.get_path('purelib'))")

WORKDIR /workdir
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/immunopipe/docker/entry.sh" ]
