# RecMpox 0.0.5 - recombination detection for mpox consensus genomes
FROM condaforge/mambaforge:latest

RUN mamba install -y -c conda-forge -c bioconda \
    python>=3.9 \
    minimap2 \
    samtools \
    squirrel \
    "iqtree=2.4.0" \
    r-base \
    r-ape \
    r-ggplot2 \
    r-phytools \
    bioconductor-ggtree \
    pip \
    setuptools \
    && mamba clean -afy

ENV PYTHONNOUSERSITE=1

WORKDIR /app
COPY . .
RUN pip install --no-deps --no-build-isolation .

ENTRYPOINT ["recmpox"]
CMD ["--help"]
