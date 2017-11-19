FROM r-base:latest

## Install packages from CRAN
RUN install2.r --error \
	-r "https://bioconductor.org/biocLite.R" \
	BiocGenerics \
	graph \
	graphite \
	igraph
	
COPY . /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts

CMD ["Rscript", "PoTRA_R.R"]