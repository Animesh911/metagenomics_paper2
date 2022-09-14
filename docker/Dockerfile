FROM continuumio/miniconda3
# NOTE: /resources will contain files available locally to the container. 
#  conda list --export >dependencies.yaml
#  docker build -t hunter .
#  docker images
#  docker rmi -f 236a5f9ffc64
# (Configuration file for fastqscreen)
WORKDIR /resources
COPY fastq_screen.conf .
#COPY /storage/aku048/nextflow/marine/adapters.fa .
#RUN bowtie2-build -f adapters.fa adapters

# NOTE: /workflow will be mounted and "overwritten" by NF on runtime
ENV HOME /work
WORKDIR ${HOME}

# Install and configure Conda stuff
COPY dependencies.yaml .
RUN apt-get update
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --add channels ursky
RUN conda config --add channels anaconda
RUN conda create --name env --file dependencies.yaml
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
# In order for scope docker.fixOwnership to work properly, procps needs to be installed in the container, which is not documented anywhere and causes a 3-day headache
RUN apt-get update && apt install -y procps g++ && apt-get clean
