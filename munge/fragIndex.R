#### Fragmentation layer ####

source("R/helperFunctions.R")

library(raster)
# ## load GLOB2009 from source
# glob <- raster(file.path(data.source,"/Glob009/","GLOBCOVER_L4_200901_200912_V2.3.tif"))
# ## Africa cropping
# library(rgdal)
# afr <- readOGR(dsn = file.path(data.source, "Africa"), layer = "AfricaUnion")
# glob.af <- crop(glob, afr, snap = "out")
# glob.m <- mask(glob.af, afr)
# 
# glob.f <- as.factor(glob.m)
# key.r <- read.csv(file.path(paste0(data.source,"/Glob009/","Globcover2009_Legend.csv")))
# key <- key.r[,c("Value", "Class")]
# key.bi <- data.frame(Value = key$Value, Class= 
#                        c(0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,NA))
# glob.r <- subs(glob.f, key.bi)
library(EBImage)

glob.arr <- raster(file.path(data.source, "GlobBi.tif"))
glob.arr <- as.Image(glob.arr)

## Stupid object
glob.bw <- bwlabel(glob.arr)
q <- as.matrix(glob.bw@.Data)
glob.arr <- raster(file.path(data.source, "GlobBi.tif"))

values(glob.arr) <- q
rm(glob.bw, q)

writeRaster(glob.arr, file.path(data.source, "Glob_BW.tif"), format = "GTiff")

rast <-raster(file.path(data.source, "Glob_BW.tif"))

Frag.Index<- focal(rast, w=matrix(1, ncol=9, nrow=9), fun=function(x){length(unique(x))} )
writeRaster(Frag.Index, file.path(data.source, "fragIndexSourceRes.tif"), format = "GTiff")

## Looks corect engough deleting progress items and crunching to study resolution in munge/geoSmash
