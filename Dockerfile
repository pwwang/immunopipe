FROM mambaorg/micromamba:latest

RUN mkdir -p /home/$MAMBA_USER/immunopipe
RUN chown $MAMBA_USER:$MAMBA_USER /home/$MAMBA_USER/immunopipe

USER $MAMBA_USER
WORKDIR /home/$MAMBA_USER/immunopipe
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

RUN micromamba install -n base --yes --file env.lock && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Install python dependencies
RUN pip install -U poetry && \
    poetry config virtualenvs.create false && \
    poetry install -v

CMD ["immunopipe"]
