# TAG dnastack/signal:e6cae1e
FROM continuumio/miniconda3:4.9.2

# Based on https://github.com/jaleezyy/covid-19-signal/blob/master/resources/Dockerfile_pipeline
MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq --allow-releaseinfo-change update && \
	apt-get -qq install \
		curl \
		wget \
		git \
		build-essential \
		libtbb2 \
		tabix \
		bcftools \
		zip \
		xxd && \
	conda init bash

# Download required datasets; these are included as part of the Dockerfile to ensure AWS's EBS does not get overwhelmed at large sample numbers
# N.B. that this makes the final Docker image ~10 GB
RUN mkdir /data
RUN	wget https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/Kraken2.tar.gz && \
	tar -zxvf Kraken2.tar.gz --directory /data && \
	rm Kraken2.tar.gz
RUN wget https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/composite_human_viral_reference.tar.gz && \
	tar -zxvf composite_human_viral_reference.tar.gz --directory /data && \
	rm composite_human_viral_reference.tar.gz

ENV KRAKEN2_DB /data/Kraken2/db
ENV COMPOSITE_REFERENCE /data/composite_human_viral_reference/composite_human_viral_reference.fna

# install snakemake and system dependencies
RUN conda create --name snakemake --channel conda-forge --channel bioconda snakemake==6.0.5 pandas

# clone the workflow repo (tag later to specific release to match software version)
RUN git clone https://github.com/jaleezyy/covid-19-signal

# set workdir to pipeline folder in repo
WORKDIR /covid-19-signal

ENV SIGNAL_VERSION e6cae1eaa375409f368aa50c6616e8bb24e0243a
RUN git checkout ${SIGNAL_VERSION}

# install all snakemake rule conda envs using dummy snakemake file
COPY ./dependencies resources/dependencies
RUN conda run \
	-n snakemake snakemake \
	--verbose \
	--use-conda \
	--conda-prefix=$HOME/.snakemake \
	-s resources/dependencies \
	--cores 8

COPY ./scripts/* /usr/local/bin/
