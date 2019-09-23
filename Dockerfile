FROM ubuntu:16.04

MAINTAINER Kristian Peters (kpeters@ipb-halle.de)

# Add cran R backport
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install apt-transport-https perl
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Generate locales
ENV LC_ALL="en_US.UTF-8"
ENV LC_CTYPE="en_US.UTF-8"
RUN locale-gen $LC_ALL
RUN dpkg-reconfigure locales

# Install packages
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install apt-transport-https make gcc gfortran g++ libblas-dev liblapack-dev libxml++2.6-dev libexpat1-dev libxml2-dev libnetcdf-dev libssl-dev r-base r-base-dev maven texlive-latex-base texlive-latex-recommended texlive-fonts-recommended git openjdk-8-jdk-headless openjdk-8-jre-headless pkg-config parallel wget curl git unzip zip python3 openbabel

# Install R packages
RUN R -e 'install.packages(c("devtools","readxl","webchem","jsonlite","rcdk","circlize","plotrix","squash"))'

# Install  Bioconductor 
RUN R CMD javareconf
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("multtest","MSnbase","mzR","MassSpecWavelet","S4Vectors","BiocStyle","faahKO","msdata","xcms","CAMERA"), ask=FALSE)'
RUN R -e 'library(devtools); devtools::install_github("https://github.com/MassBank/RMassBank", ref="treutler-merge")'
RUN R -e 'library(devtools); install_github(repo = "CDK-R/rinchi@master")'

# Fix RMassBank
WORKDIR /usr/src
RUN git clone --branch treutler-merge https://github.com/MassBank/RMassBank

# Cleanup
RUN apt-get -y --purge --auto-remove remove make gcc gfortran g++
RUN apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*

# Add scripts folder to container
ADD galaxy/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*

