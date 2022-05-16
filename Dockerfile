FROM continuumio/miniconda3:4.11.0

WORKDIR /immunopipe
COPY . /immunopipe

# Install dependencies
RUN conda env create --file docker/environment.yml && \
    conda clean --all --yes

# See https://github.com/python-poetry/poetry/issues/3628
# for the patching poetry (that's why we have to pin poetry's version)
# The python dependencies are installed by conda,
# use poetry to install again to check they are at the right versionsls /ho
SHELL [ "/bin/bash", "--login", "-c" ]
RUN conda activate immunopipe && \
    python -m pip install poetry==1.1.13 && \
    python -m poetry config virtualenvs.create false && \
    patch $(python -c 'from poetry.repositories import installed_repository as ir; print(ir.__file__, end="")') docker/poetry.pth && \
    python -m poetry install -v && \
    # For singularity to init conda
    mkdir /workdir && \
    mkdir -p /home/immunopipe_user/ && \
    cp docker/.bashrc /home/immunopipe_user/

ENTRYPOINT [ "conda", "run", "-n", "immunopipe", "python", "-m", "immunopipe" ]
