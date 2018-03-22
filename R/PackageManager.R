#### Package Manager ####
#### FS Sept 2017    ####
#########################

require(devtools)
local({r <- getOption("repos")
r["CRAN"] <- "https://cran.r-project.org"
options(repos=r)
})
getOption("repos")

packs <- rownames(installed.packages()) # names of all installed packages
req <- c("biomod2", "raster", "rgeos", "rgdal", "RStoolbox",
         "dplyr", "tidyr", "snowfall", "snow", "gtools", "spatstat",
         "skimr", "ggplot2", "rlang", "lazyeval", "gridExtra")

toInstall <- req[!is.element(req,packs)] # packages needing installation

if(length(toInstall)>0){
  install.packages(toInstall)
}

# clean up workspace:
rm(packs,req,toInstall)
