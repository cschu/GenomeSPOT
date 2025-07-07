FROM ubuntu:22.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="0.1"
LABEL description="This is a Docker image for GenomeSPOT"


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt-get install -y wget python3-pip python-is-python3 git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev rsync

RUN apt clean

ADD LICENSE README.md requirements.txt setup.cfg setup.py /opt/software/genomespot/
ADD genome_spot /opt/software/genomespot/genome_spot/
ADD data /opt/software/genomespot/data/
ADD models /opt/software/genomespot/models/

RUN cd /opt/software/genomespot/ && \
	pip install . && \
	pip install -r requirements.txt
