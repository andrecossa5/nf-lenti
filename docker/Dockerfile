# MI_TO image

################################################################

## Start Debian basic image: NBBB ricordarsi --platform linux/x86_64 in docker build
FROM python:3.9-bullseye

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
    && apt-get --quiet upgrade --yes \
    && apt-get --quiet install --yes \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      libz-dev \
      wget \
      procps \
      git \
      cmake \
    && rm -rf /var/lib/apt/lists/*

## Mamba
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    && bash Mambaforge-Linux-x86_64.sh -b -p /mambaforge
ENV PATH="/mambaforge/bin/:$PATH"

# Tools
RUN mamba install -c bioconda samtools cellsnp-lite freebayes
RUN git clone https://github.com/caleblareau/maegatk.git
RUN git clone https://github.com/KChen-lab/Monopogen.git \
   && cd Monopogen \
   && pip install -e .


## Default command to initialize the container
CMD [ "bash" ]

#################################################################

