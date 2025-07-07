FROM ubuntu:22.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="0.1"
LABEL description="This is a Docker image for GenomeSPOT"


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt-get install -y wget python3-pip python-is-python3 git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev rsync

RUN apt clean

RUN mkdir -p /opt/software
WORKDIR /opt/software

RUN git clone https://github.com/cultivarium/GenomeSPOT.git && \
    cd GenomeSPOT && \
    pip install . && \
    pip install -r requirements.txt

