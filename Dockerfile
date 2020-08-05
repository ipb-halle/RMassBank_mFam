FROM ubuntu:20.04

MAINTAINER Kristian Peters (kpeters@ipb-halle.de)

# Add cran R backport
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install ca-certificates apt-transport-https perl locales gnupg2
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# Generate locales
ENV LC_ALL="en_US.UTF-8"
ENV LC_CTYPE="en_US.UTF-8"
RUN locale-gen $LC_ALL
RUN dpkg-reconfigure locales

# Install packages
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y --allow-unauthenticated install apt-transport-https make gcc gfortran g++ libblas-dev liblapack-dev libxml++2.6-dev libexpat1-dev libxml2-dev libnetcdf-dev libssl-dev pkg-config wget curl git unzip zip python3 python3-pip r-base r-base-dev pandoc pandoc-data openjdk-11-jdk pkg-config parallel wget curl git unzip zip python3 openbabel maven

# Install R packages
RUN R -e 'install.packages(c("devtools","readxl","webchem","jsonlite","rcdk","circlize","plotrix","squash"))'

# Install  Bioconductor 
RUN R CMD javareconf
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(ask = FALSE, update = TRUE); BiocManager::install(c("multtest","MSnbase","mzR","MassSpecWavelet","S4Vectors","BiocStyle","faahKO","msdata","xcms","CAMERA"), ask = FALSE, update = TRUE)'
RUN R -e 'library(devtools); devtools::install_github("https://github.com/MassBank/RMassBank", ref="treutler-merge")'
#RUN R -e 'library(devtools); devtools::install_github("https://github.com/MassBank/RMassBank")'
RUN R -e 'library(devtools); install_github(repo = "CDK-R/rinchi@master")'

# Install old version of webchem
RUN R -e 'devtools::install_github("ropensci/webchem", ref="8c9324225d161158bf12e26f688a331201ec9647")'

# Install RMassBank
WORKDIR /usr/src
RUN git clone --branch treutler-merge https://github.com/MassBank/RMassBank
#RUN git clone https://github.com/MassBank/RMassBank

# Install MassBank
WORKDIR /usr/src
RUN git clone https://github.com/MassBank/MassBank-web
WORKDIR /usr/src/MassBank-web/MassBank-Project/
RUN mvn package

# Cleanup
RUN apt-get -y --purge --auto-remove remove make gcc gfortran g++
RUN apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*

# Add scripts folder to container
ADD galaxy/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*

