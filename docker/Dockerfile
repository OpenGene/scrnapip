FROM rocker/rstudio:latest
 
MAINTAINER Li yq <liyq@haplox.com>
 
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
 libxtst-dev libpng++-dev libxml2-dev cmake libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libgit2-dev\
  && apt-get clean all && \
  apt-get purge &&\
  mkdir /home/dataa
#COPY ./data/icu4c-69_1-data-bin-l.zip /home/dataa/icu4c-69_1-data-bin-l.zip

RUN Rscript -e "options(timeout = 1500);options(repos=structure(c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')));install.packages(c('rmarkdown', 'BiocManager','devtools'));install.packages('stringi');options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');install.packages('configr');install.packages('Seurat')"

RUN Rscript -e "options(timeout = 1500);options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');options(repos=structure(c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/'))); BiocManager::install(c('SingleR','monocle','biomaRt','ggtree','GSVA','qvalue'));install.packages('optparse');install.packages('logging');install.packages('circlize')"


RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y && \
  apt-get install -y libcairo2-dev libxt-dev libnetcdf-dev python3-pip  &&\ 
  apt-get install -y libgdal-dev libproj-dev libgeos-dev libudunits2-dev libnode-dev libcairo2-dev libnetcdf-dev &&\
  apt-get install -y libglu1-mesa-dev freeglut3-dev mesa-common-dev &&\
  ln -s /usr/bin/python3 /usr/bin/python &&\
  pip3 install numpy -i https://pypi.tuna.tsinghua.edu.cn/simple &&\
  pip3 install ScanoramaCT -i https://pypi.tuna.tsinghua.edu.cn/simple 
RUN Rscript -e "options(timeout = 1500);options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');options(repos=structure(c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')));install.packages(c( 'colourpicker', 'formattable', 'msigdbr','readr', 'shinydashboard','shinyWidgets','DT'));BiocManager::install('clusterProfiler');BiocManager::install(c('pathview','org.Hs.eg.db','org.Mm.eg.db','Cairo','topGO'));install.packages(c('harmony','R.utils','checkmate','pROC'));BiocManager::install('scran');BiocManager::install('ExperimentHub'); BiocManager::install('genomicInstability'); BiocManager::install('ComplexHeatmap'); BiocManager::install('sva')"

COPY ./data /home/dataa

RUN Rscript -e "options(timeout = 1500);options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');options(repos=structure(c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')));BiocManager::install('batchelor');install.packages('/home/dataa/cerebroApp-1.2.2.tar.gz',type='source');install.packages('/home/dataa/monocle/monocle_2.28.0.tar.gz',type='source');options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');options(repos=structure(c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/'))); install.packages('/home/dataa/dlm_1.1-5.tar.gz',repos = NULL);install.packages('/home/dataa/reactome.db_1.82.0.tar.gz',repos = NULL) ;BiocManager::install('ReactomePA');devtools::install_local('/home/dataa/seurat-wrappers-master.zip',upgrade ='never'); devtools::install_local('/home/dataa/copykat-master.zip',quiet =T,upgrade ='never');devtools::install_local('/home/dataa/CytoTRACE_0.3.3.tar.gz',upgrade ='never');devtools::install_github('sqjin/CellChat') "

RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
