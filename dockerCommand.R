devtools::install_github("o2r-project/containerit")
install.packages("rgdal")
install.packages("maptools")
library("containerit")
library("rgdal")
library("maptools")

# create simple script file
scriptfile <- tempfile(pattern = "containerit_", fileext = ".R")
writeLines(c('library(rgdal)',
             'nc <- rgdal::readOGR(system.file("shapes/", package="maptools"), "sids", verbose = FALSE)',
             'proj4string(nc) <- CRS("+proj=longlat + datum=NAD27")',
             'plot(nc)'), scriptfile)

# use a custom startup command
scriptCmd <- CMD_Rscript(basename(scriptfile))

# create Dockerfile for the script
dockerfile_object <- dockerfile(from = scriptfile, silent = TRUE, cmd = scriptCmd)

print(dockerfile_object)