##Geo munge from source 
source("R/helperFunctions.R")
## African extent
library(rgdal); library(rgeos); library(maptools)
# afr <- readOGR(dsn = path.expand(file.path(data.source, "Africa")),
#                layer = "AfricanCountires") 
# afr <- gUnaryUnion(afr, id = afr@data$CONTINENT)
# afr.sp <- as(afr, "SpatialPolygonsDataFrame")
# class(afr.sp)
# 
# writeOGR(afr.sp, dsn = path.expand(file.path(data.source, "Africa")),
#          layer = "AfricaUnion", driver = "ESRI Shapefile")

# afr.u <- readOGR(dsn = file.path(data.source, "Africa"),
#                  layer = "AfricaUnion")

##
library(raster)

  ##Method idea from Fhar 2016
# prec.wc <- file.path(data.source, "wc2.0_5m")
# prec.ls <- lapply(paste0(prec.wc, "/", list.files(prec.wc)), raster)
# prec.stk <- do.call(stack, prec.ls)
# prec.ann <- sum(prec.stk)
# #Create mask
# m <-c(-1,500,NA, 500,cellStats(prec.ann, max),1)
# rf.mask.500 <- reclassify(prec.ann,m)
# rf.25 <- projectRaster(rf.mask.500, res = c(res(rf.mask.500)*2.5),
#                        crs = proj4string(rf.mask.500), method = "ngb")
# rf.c <- crop(rf.25, afr.u, snap = "out")
# rf <- mask(rf.c, afr.u)
# writeRaster(rf, paste0(data.source,"/cropMask.tif"), overwrite = T)
crop.mask <- raster(paste0(data.source,"/cropMask.tif"))

## RainFall Data Sourced from worldclim
prec.wc <- file.path(data.source, "wc2.0_5m")
prec.ls <- lapply(paste0(prec.wc, "/", list.files(prec.wc)), raster)
prec.stk <- do.call(stack, prec.ls)
names(prec.stk) <- paste0("prec",1:12)
prec.clean <- mask(projectRaster(prec.stk, crop.mask), crop.mask)
writeRaster(prec.clean, paste0(clean.dir,"/", names(prec.clean), ".tif"), 
            bylayer = T, overwrite = T)
rm(prec.wc, prec.ls, prec.stk, prec.clean)

## Mean Temp Data sourced from world clim
tmean <- file.path(data.source, "tmean")
tmean.ls <- lapply(paste0(tmean,"/",list.files(tmean, pattern = "*.bil")),raster)
tmean.stk <- do.call(stack, tmean.ls)
tmean.clean <- mask(projectRaster(tmean.stk, crop.mask), crop.mask)
writeRaster(tmean.clean, paste0(clean.dir,"/", names(tmean.clean),".tif"),
            bylayer = T, overwrite = T)
rm(tmean, tmean.ls, tmean.stk, tmean.clean)

## Bclim data sourced from world clim
bclim <- file.path(data.source, "bio")
bio.ls <- lapply(paste0(bclim,"/",list.files(bclim, pattern = "*.bil")),raster)
bio.stk <- do.call(stack, bio.ls)
bio.clean <- mask(projectRaster(bio.stk, crop.mask),crop.mask)
writeRaster(bio.clean, paste0(clean.dir,"/", names(bio.clean),".tif"),
            bylayer = T, overwrite = T)
rm(bclim, bio.ls, bio.stk, bio.clean)

## Diversity rasters Data generated from the IUCN Data set (see gen script)
div <- file.path(data.source, "DiversityRasters")
div.ls <- lapply(paste0(div,"/",list.files(div)),raster)
div.stk <- do.call(stack, div.ls)
div.clean <- mask(projectRaster(div.stk, crop.mask), crop.mask)
writeRaster(div.clean, paste0(clean.dir,"/", names(div.clean),".tif"),
            bylayer = T, overwrite = T)
rm(div, div.ls, div.stk, div.clean)

## Potential Evapotransporation See read me for info
pet <- file.path(data.source, "PET_he_monthly")
pet.ls <- list()
for(i in 1:12){
  tp <- raster(paste0(pet, "/pet_he_",i,"/w001001.adf"))
  names(tp) <- paste0("pet",i)
  dnc <- mask(projectRaster(tp, crop.mask),crop.mask)
  writeRaster(dnc, paste0(clean.dir,"/", names(dnc),".tif"), overwrite = T)
}
rm(pet, pet.ls, i, dnc, tp)

## Enhanced Vegatitive Index From Modis (see read me)
evi <- file.path(data.source,"EVI")
for(i in 1:12){
  files <- c(paste0(i,"eviSA.tiff"),paste0(i,"eviWA.tiff"),paste0(i,"eviEA.tiff"))
  evi.l <- lapply(paste0(evi,"/",files), raster)
  evi.a <- do.call(merge, evi.l)
  names(evi.a)  <- paste0("evi",i)
  dnc <- mask(projectRaster(evi.a, crop.mask),crop.mask)
  writeRaster(dnc, paste0(clean.dir,"/", names(dnc),".tif"), overwrite = T)
}

rm(files,evi.l,evi.a,dnc)

## LandCover Data sourced from GlobCover2009
glob <- raster(paste0(data.source,"/Glob009/","GLOBCOVER_L4_200901_200912_V2.3.tif"))
glob.f <- as.factor(glob)
key.r <- read.csv(file.path(paste0(data.source,"/Glob009/","Globcover2009_Legend.csv")))
key <- key.r[,c("Value", "Class")]
glob.m <- projectRaster(glob.f, crop.mask, method = "ngb")
glob.mm <- mask(glob.m, crop.mask)
glob.r <- subs(glob.mm, key)
writeRaster(glob.r, file.path(clean.dir,"LandCover.tif"),format = "GTiff", overwrite = T)

## population Density 
pop.den <- raster(file.path(data.source, "af_gpwv3_pdens_00_ascii_25", "afds00ag.asc"))
proj4string(pop.den) <- proj4string(crop.mask)
pop.prj <- projectRaster(pop.den, crop.mask)
pop.m <- mask(pop.prj, crop.mask)
writeRaster(pop.m, file.path(clean.dir,"popDen.tif"), format = "GTiff")

## Fragmentation index
frag <- raster(file.path(data.source, "fragIndexSourceRes.tif"))
frag.prj <- projectRaster(frag, crop.mask)
frag.m <- mask(frag.prj, crop.mask)
writeRaster(frag.m, file.path(clean.dir, "fragIndex.tif"), format = "GTiff")
