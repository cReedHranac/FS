#### Animal Oubreak point fix ####

rm(list = ls())
source("R/helperFunctions.R")
library(data.table); library(dplyr)

dat <- tbl_df(fread(file.path(clean.dir, "longTable.csv")))
prob <- dat %>% 
  filter(cell == "c65880" )
prob$x
prob$y
## Mystery occurence in this position. #WTF
#lets see if it falls in the cropping bounds we use. 
library(raster)
rf.mask <- raster(file.path(data.source, "cropMask.tif"))
cellFromXY(rf.mask, c(prob$x[[1]], prob$y[[1]])) # same
rf.mask[cellFromXY(rf.mask, c(prob$x[[1]], prob$y[[1]]))] 
#outside of cropping area... But should it be?

### Extracted from ViralSmash
library(rgdal);library(rgeos);library(raster);library(dplyr)

an.src <- read.csv(file.path(data.source, "vIndex", "latest_animal_case.csv"))
an.simple <- an.src %>% dplyr::select(OUTBREAK_ID, Month.Start)
colnames(an.simple)[1] <- "Outbreak_ID"

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

ma.an <- inner_join(as.data.frame(an.ma), an.simple, by = "Outbreak_ID")

coordinates(ma.an) <- ~ x+y
plot(rf.mask)
points(ma.an)
points(xyFromCell(rf.mask,65880), col = "red")

### Stupid freaking thing is in the ocean.... are there others out there too?
rf.mask[cellFromXY(rf.mask,ma.an@coords)] ## appears to only be that one...
rowColFromCell(rf.mask, 65880)# going to try one cell north
rf.mask[cellFromRowCol(rf.mask, 154, 122)] # that appears to work. 
new.cell <- cellFromRowCol(rf.mask, 154, 122)
points(xyFromCell(rf.mask,new.cell), col = "blue")
## all appears correct

new.xy <- xyFromCell(rf.mask, 65453) #extract new xy
## Inserted into dataframe via Viral Smash, re-source ViralSmash and FrameGen
