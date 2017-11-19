FROM r-base:latest

# install R package dependencies
RUN apt-get update && apt-get install -y \
    bioconductor/release_base2:latest
	
COPY . /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts

CMD ["Rscript", "PoTRA_R.R"]