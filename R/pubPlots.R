###############################################################################
#### Figures for publication  29/01/2018                               ########
###############################################################################
source("R/helperFunctions.R")

library(ggplot2); library(dplyr); library(data.table)
#### Figure 1 ####
#3 pannel figure, 1 pannel, map w/ outbreak locations (symbols for an/hum), 
# and 2 outbreak time lines, 1 entire history, one collapsed year
  #### data ####
an.ob <- fread(file.path(clean.dir, "annOB.PPM.csv")) #animal points
an.an <- fread(file.path(data.source, "Animal_Index_Jan2018.csv")) #notes on outbreaks
an.sub <- an.an %>%
  dplyr::select(OUTBREAK_ID, Org_smp, Year_Start)
names(an.sub)[1] <- names(an.ob)[1]
an.full <- an.sub %>%
  left_join(an.ob,an.sub, by = "Outbreak_ID")

hum.ob <- fread(file.path(clean.dir, "humOB.PPM.csv"))
hum.an <- fread(file.path(data.source, "Human_Index_30May.csv"))
hum.sub <- hum.an %>%
  dplyr::select(Outbreak_ID, Year_Start)
hum.full <- hum.sub %>%
  left_join(hum.ob, hum.sub, by = "Outbreak_ID")


hum.full$Org_smp <- "human"
hum.s <- hum.full %>%
  dplyr::select(which(names(hum.full) %in% names(an.full)))

ob.full <- bind_rows(hum.s, an.full)
colnames(ob.full)[3:4] <- c( "long", "lat")
ob.full$Org_smp <- as.factor(ob.full$Org_smp)

z <- list()
for(i in 1:nrow(ob.full)){
  ifelse(ob.full[i,"Org_smp"] == "human", z[[i]] <- "human",z[[i]] < "animal")  
}

# library(devtools)
# devtools::install_github("sckott/rphylopic")
library(rphylopic)
## name id ##
bat <- 	104257
chimpanzee <- 109082
duiker <-  2478895
gorilla <-  763767
human <- 	109086

critters <- c(bat, chimpanzee, gorilla, duiker, human)
critters.id <- lapply(critters, ubio_get)
critter.names <- lapply(critters.id, name_images) ## Needs manual attention
c.names <- c(critter.names[[1]]$same[[1]]$uid, 
             critter.names[[2]]$same[[1]]$uid,
             critter.names[[3]]$same[[1]]$uid,
             critter.names[[4]]$supertaxa[[1]]$uid,
             critter.names[[5]]$same[[1]]$uid)
critter.img <- lapply(c.names, image_data, size = 64)
c.img <- list()
for(i in 1:length(c.names)){
  c.img[[i]] <- critter.img[[i]][[1]]
}

## create vector of images
img.list <- list()
for(i in 1:nrow(ob.full)){
  j <- match(ob.full[i,"Org_smp"], levels(ob.full$Org_smp)) 
  img.list[[i]] <- c.img[[j]]
}

## create color vector
cz <- c("darkorange2", "black", "green3", "dodgerblue2", "magenta1")
col.list <- list()
for(i in 1:nrow(ob.full)){
  j <- match(ob.full[i,"Org_smp"], levels(ob.full$Org_smp)) 
  col.list[[i]] <- cz[[j]]
}

#### Pannel 1 Map ####
library(rgdal);library(rgeos);library(raster)
## accessory layers
afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")
rf.poly <- rasterToPolygons(raster(file.path(data.source, "cropMask.tif")),
                            fun = function(x){x==1}, dissolve = T)
bkg <- theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  plot.title = element_text(hjust = 0.5))

ob.plot <- ggplot() +
  geom_polygon(data = fortify(afr.poly),
               aes(long, lat, group = group), 
               colour = "grey20",
               alpha = .25) +
  coord_fixed()+
  aes(x=long, y=lat) +
  geom_polygon(data = fortify(rf.poly),
               aes(long, lat, group = group),
               colour = "white", 
               alpha = .25,
               fill = "cornsilk") +
  coord_cartesian(xlim = c(-20, 53),
                  ylim = c(-36, 40))
  
  
for( i in 1:nrow(ob.full)){
  ob.plot <- ob.plot + add_phylopic(img.list[[i]], 1, 
                                    ob.full$long[i],
                                    ob.full$lat[i],
                                    ysize = 2.5,
                                    color = col.list[[i]])
}
  
  
ob.plot + bkg
# geom_point(data = fortify(ob.full),
#            aes(x = long, y = lat,
#                group = Org_smp,
#                shape = Org_smp,
#                colour = Org_smp)) +
# scale_shape_manual(values = c(15:19))+
# add_phylopic(c.img[[3]])
