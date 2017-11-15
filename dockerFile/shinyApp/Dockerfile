FROM rocker/shiny
MAINTAINER Wentao Zhou

# install R package dependencies
RUN apt-get update && apt-get install -y \
    libssl-dev \
    ## clean up
    && apt-get clean \ 
    && rm -rf /var/lib/apt/lists/ \ 
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
	
## Install packages from CRAN
RUN install2.r --error \ 
    -r 'http://cran.rstudio.com' \
	googleAuthR \
	devtools \
    && Rscript -e "devtools::install_github(c('wtzhou74/dockerPotra'))" \
    ## clean up
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

## assume shiny app is in build folder /shiny
COPY ./shinyApp/ /srv/shiny-server/shinyapp/