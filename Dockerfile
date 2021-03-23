### copyright 2017-2021 Regents of the University of California and the Broad Institute. All rights reserved.

FROM genepattern/seurat-suite:2.4
# The Dockerfile for this is on the folder named Docker
MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C
USER $NB_USER

RUN R -e "library('Seurat');sessionInfo()"

#RUN mkdir /module
#RUN mkdir /temp
ADD src /module/
RUN ls /module/

# build using this:
# docker build -t genepattern/seurat-batchcorrection:3.1 .
