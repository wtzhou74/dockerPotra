FROM rocker/shiny
MAINTAINER Wentao Zhou
    
## Install packages from CRAN
RUN install2.r --error \ 
    -r 'http://cran.rstudio.com' \
	googleAuthR \
    && Rscript -e "devtools::install_github(c('wtzhou74/dockerPotra')" \
    ## clean up
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

## assume shiny app is in build folder /shiny
COPY ./shinyApp/ /srv/shiny-server/shinyapp/