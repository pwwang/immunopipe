FROM continuumio/miniconda3:4.11.0

WORKDIR /immunopipe
COPY . .

# Mostly R and R packages
# .biopipen.toml specifies the Rscript in the environment
RUN conda env create --file environment.yml && \
    conda clean --all --yes

# See https://github.com/python-poetry/poetry/issues/3628
# for the patching
SHELL [ "/bin/bash", "--login", "-c" ]
RUN conda activate immunopipe && \
    python -m pip install poetry==1.1.13 && \
    python -m poetry config virtualenvs.create false && \
    patch $(python -c 'from poetry.repositories import installed_repository as ir; print(ir.__file__, end="")') poetry.pth && \
    python -m poetry install -v

ENTRYPOINT [ "conda", "run", "-n", "immunopipe", "python", "-m", "immunopipe" ]
