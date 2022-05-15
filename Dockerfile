FROM continuumio/miniconda3:4.11.0

WORKDIR /immunopipe
COPY . .

# Mostly R and R packages
# .biopipen.toml specifies the Rscript in the environment
RUN conda env create --file environment.yml && \
    conda clean --all --yes

SHELL [ "/bin/bash", "--login", "-c" ]
RUN conda activate immunopipe && \
    python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v

ENTRYPOINT [ "conda", "run", "-n", "immunopipe", "python", "-m", "immunopipe" ]
