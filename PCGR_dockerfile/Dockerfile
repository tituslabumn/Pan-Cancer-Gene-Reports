FROM rocker/r-base:latest

RUN mkdir PCGR

COPY startup.sh STARTUP.txt dockerfile_Rpackages.R /

RUN apt-get update && apt-get install -y \
nano \
libxml2-dev \
libcurl4-openssl-dev \
libssl-dev \
pandoc \
default-jdk \
git 

RUN Rscript dockerfile_Rpackages.R

ENTRYPOINT /./startup.sh && /bin/bash
WORKDIR /PCGR/Pan-Cancer-Gene-Reports/

