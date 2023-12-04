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

## Clone repo and install code
RUN echo new_commits
RUN git clone https://github.com/andrecossa5/mito_utils.git \
    && cd mito_utils \
    && mamba env update -n base -f envs/environment.yml \
    && mamba develop .. \
    && mamba develop ../mito_utils/ 

# pysam dnaio
RUN pip install pysam dnaio

# Maegatk repo
RUN git clone https://github.com/caleblareau/maegatk.git

## Default command to initialize the container
CMD [ "bash" ]

#################################################################

