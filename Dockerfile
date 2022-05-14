FROM continuumio/conda-ci-linux-64-python3.9:latest

WORKDIR /immunopipe
COPY . /immunopipe/

RUN conda create -n mamba -c conda-forge mamba
RUN conda activate mamba
RUN mamba create -n immunopipe -f environment.yml
RUN conda activate immunopipe
RUN pip install -U poetry
RUN poetry config virtualenv.create false
RUN poetry install -v

ENTRYPOINT ["./entrypoint.sh"]
