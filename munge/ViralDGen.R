#### Viral Occurrence dataframes Munge ####
#### April 2018                        ####
###########################################

## script for creating the dataframes to be included as SI materials with publication
## Largely taken and modified from munge/ViralSmash.R

source("R/helperFunctions.R")
library(rgdal);library(rgeos);library(raster);library(dplyr)

c.mask <- raster(file.path(data.source, "cropMask.tif"))
wgs <- proj4string(c.mask)

#### human index cases ####
## note : bring in the Human_Index_30May. Contains 1 additional point over
## the source dataframe provided representing the Inkanamo outbreak
## UPDATE 12.2.018: additional entry for 2017 DRC outbreak

hum.src <- read.csv(file.path(data.source, "vIndex", "Human_Index_12_2_18.csv"))

#### Polygons ####
hum.poly <- readOGR(file.path(data.source, "vIndex", "humanORG"), 
                    "index_human_case_modified_polygon")
##addition of outbreak 32
drc.poly <- readOGR(file.path(data.source, "zone_ste_puc"), 
                    "Zone_Ste_Puc")
##pull out the zone sante
likati <- drc.poly[drc.poly$Nom_ZS_PUC == "Likati",]
lik.prj <- spTransform(likati, proj4string(hum.poly))
lik.prj$OUTBREAK <- 32

##bind
hum.poly <- rbind(hum.poly, lik.prj[,"OUTBREAK"])


## Extract centroids but maintain the data from the original obj
hum.poly.c <- SpatialPointsDataFrame(gCentroid(hum.poly, byid = T),
                                     hum.poly@data, match.ID = F)
#### Points ####
## Still needs to be added here
hum.pts <- readOGR(file.path(data.source, "vIndex", "humanORG"), 
                   "index_human_case_modified_point")
new <- c(-0.616667, 20.433333, 31, 20.433333, -0.616667) ## outbreak 31
human.add <- rbind(as.data.frame(hum.pts), new)
coordinates(human.add) <- ~ coords.x1 + coords.x2
proj4string(human.add) <- CRS(wgs)
## Remove the lat and long cols
hum.add <- human.add[,3]

#### masterFrame ####
hu.ma <- rbind(hum.poly.c, hum.add)
names(hu.ma) <- "Outbreak_ID"

## Combine the exracted poitns with the source dataframe 
ma.hu <- inner_join(as.data.frame(hu.ma), hum.src, by = "Outbreak_ID")

## Clean them 
human.out <- ma.hu[,colnames(ma.hu) %!in% c("Lat", "Lon")]
colnames(human.out)[2:3] <- c("lon", "lat")

#write out
write.csv(human.out, "data/HumOutbreakDB.csv", row.names = F)


#### Animal Index cases ####
an.src <- read.csv(file.path(data.source, "vIndex", "Animal_Index_Jan2018.csv"))
colnames(an.src)[2] <- "Outbreak_ID"

## remove the sero data that does not fit the analysis (ID 45:50)
an.src <- an.src[an.src$Outbreak_ID %!in% 45:50,]

#### Polygons ####
an.poly <- readOGR(file.path(data.source, "vIndex", "animalORG"),
                   "animal_polygon")
an.poly.c <- SpatialPointsDataFrame(gCentroid(an.poly, byid = T),
                                    an.poly@data, match.ID = F)
names(an.poly.c) <- "Outbreak_ID"

#### Points ####
an.pt <- readOGR(file.path(data.source, "vIndex", "animalORG"),
                 "animal_points")
an.pt <- an.pt[,1]
names(an.pt) <- "Outbreak_ID"

#### masterFrame ####
an.ma <- rbind(an.poly.c, an.pt)
ma.an <- inner_join(as.data.frame(an.ma), an.src, by = "Outbreak_ID")

## clean for write out
not.needed <- c("Lat.2", "Long.2", "Lat.3", "Long.3",
                "Lat.4", "Long.4", "Lat.5", "Long.5",
                "Lat.6", "Long.6", "Lat.7", "Long.7",
                "Lat.8", "Long.8", "Lat.9", "Long.9",
                "Lon", "Lat", "Map.to.digitize","Checked", "X", "old_id")
an.out <- ma.an[,colnames(ma.an) %!in% not.needed]

write.csv(an.out, "data/AnOutbreakDB.csv", row.names = F)
