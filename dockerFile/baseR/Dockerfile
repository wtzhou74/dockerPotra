FROM r-base:latest

RUN Rscript -e 'source("https://bioconductor.org/biocLite.R");biocLite(c("BiocStyle","graph","graphite","igraph"))'
	
COPY . /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts

CMD ["Rscript", "PoTRA_R.R"]