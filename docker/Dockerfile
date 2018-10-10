FROM r-base:3.1.3

COPY sources.list /etc/apt/sources.list

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install build-essential --yes && \
    apt-get install python-dev --yes && \
    apt-get install default-jre --yes && \
    wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py 

RUN pip install numpy

#
# install R packages 
#
#   from the Broad dotkits circa DAPPLE v0.17 
#     ssh gpprod@gpprod
#     use Python-2.6
#     use R-3.1
#     R
#     > library(plotrix)
#     > library(igraph)
#     > sessionInfo()
#     R version 3.1.1 (2014-07-10)
#     Platform: x86_64-unknown-linux-gnu (64-bit)
#     
#     locale:
#     [1] C
#     
#     attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
#     
#     other attached packages:
#     [1] igraph_0.7.1   plotrix_3.5-12
#

# add this line to the Rprofile.site, needed to install packages from https sites
#   options(download.file.method = "wget")
# 
RUN echo "\\noptions(download.file.method = \"wget\")\\n" >> /etc/R/Rprofile.site

RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/plotrix/plotrix_3.5-12.tar.gz', \
  repos = NULL, type='source')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz', \
  repos = NULL, type='source')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_0.7.1.tar.gz', \
  repos = NULL, type='source')"

# copy data into container
COPY source/code /code
COPY source/data /data

# run shell script to download the data
RUN /code/downloadData.sh
