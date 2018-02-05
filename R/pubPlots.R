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

#### Pannel 1 Maps ####
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
  geom_polygon(data = fortify(rf.poly),
               aes(long, lat, group = group),
               colour = "white", 
               alpha = .25,
               fill = "cornsilk")+
  coord_fixed(xlim = c(-20, 53),ylim = c(-36, 40))
  
  
for( i in 1:nrow(ob.full)){
  if(ob.full$Org_smp[i] == "human"){
    j <- 3
  } else{ifelse(ob.full$Org_smp[i] == "bat",  j <- 1.9,  j <- 2.5)}
  ob.plot <- ob.plot + 
    add_phylopic(img.list[[i]], .5, 
                 ob.full$long[i],
                 ob.full$lat[i],
                 ysize = j,
                 color = col.list[[i]]) 
    
}

ob.plot + bkg ## good enough for now...


### Insert plot
## Need for rain mask? just shows low resolution
ob.insert <- ggplot()+
  geom_polygon(data = fortify(afr.poly),
               aes(long, lat, group = group), 
               colour = "white",
               alpha = .25,
               fill = "grey20") +
  coord_fixed(xlim = c(9, 17.5),ylim = c(2,-2.5 ))

for( i in 1:nrow(ob.full)){
  if(ob.full$Org_smp[i] == "human"){
    j <- .3
  } else{ifelse(ob.full$Org_smp[i] == "bat",  j <- .2,  j <- .25)}
  ob.insert <- ob.insert + 
    add_phylopic(img.list[[i]], .5, 
                 ob.full$long[i],
                 ob.full$lat[i],
                 ysize = j, ##humans need bigger, bats smaller
                 color = col.list[[i]]) 
}

ob.insert +bkg ## That's pretty good... 

#### Pannel 2 Time line ####
ob.T <- ob.full
library(zoo)
ob.T$Date <- as.yearmon(paste(ob.T$Year_Start, ob.T$Month.Start), "%Y %m")
ob.a <- ob.T %>%
  dplyr::arrange(Date)
ob.a$Month.Start <- as.factor(ob.a$Month.Start)

## Build offset size
off.list <- list()
img.T <- list()
col.T <- list()
human.T <- list()
for(i in 1:nrow(ob.T)){
  j <- match(ob.a[i,"Org_smp"], levels(ob.a$Org_smp)) 
  off.list[[i]] <- j ##Offset length
  img.T[[i]] <- c.img[[j]] ## which image to use
  col.T[[i]] <- cz[[j]] ##which color to use 
  ifelse(j == 5, human.T[[i]] <- 1, human.T[[i]] <- -1) #directionality
}
ob.a$offset <- unlist(off.list) * unlist(human.T)

## Plot 
## remove vertical lines 
g.time <- ggplot()+
  geom_segment(aes(x = Date, y = offset, xend = Date, color = Org_smp),
               data = ob.a,yend = 0)+
  geom_segment(aes(x = 1975, y = 0, xend = 2018, yend = 0),
               data = ob.a, arrow = arrow(length =  unit(x = 0.2,units = 'cm'),type = 'closed')) +
  scale_x_yearmon(format = "%Y %m", n = 10) 
# 
# for(i in 1:nrow(ob.a)){
#   g.time <- g.time +
#     add_phylopic(img.T[[i]], 1,
#                  ob.a$Date[i],
#                  ob.a$offset[i],
#                  # ysize = .3,
#                  color = col.T[[i]])
# }
g.time
#### Pannel 3 Bar with Smooth ####
### If we can start this at setptember we can likely get the bimodal distribution
### to show up better. Do two different smooths for human/nonhumans and maybe one for total?
human.B <- list()
for(i in 1:nrow(ob.a)){
  j <- match(ob.a[i,"Org_smp"], levels(ob.a$Org_smp)) 
  ifelse(j == 5, human.B[[i]] <- "Human", human.B[[i]] <- "Animal") #directionality
}
ob.a$Org <- unlist(human.B)
ob.bar <- ob.a %>% 
  group_by(Org,Month.Start) %>% 
  summarise(num=n())
ob.human <- ob.bar %>%
  filter(Org == "Human")
ob.human$Month.Start <- as.integer(ob.human$Month.Start)
ob.animal <- ob.bar %>%
  filter(Org == "Animal")
ob.animal$Month.Start <- as.integer(ob.animal$Month.Start)

## issues: how do I get time to recirculate? is it worth starting the time at 9 to
## demonstarte the bimodality of the set? Get error bars different colors?

ggplot(data=ob.a)+ geom_bar(aes(x=Month.Start,fill=Org_smp))+
  geom_smooth(data = ob.human,aes(x=Month.Start,y=num, color = Org)) +
  geom_smooth(data = ob.animal,aes(x=Month.Start,y=num, color = Org))+
  bkg



#### Figure 2 ####
# Map pannel(s) with bat locations, pannel with collaped year & smoots for 
# each class like above (may need to be seperated by classes if too busy)

  ### Data #### 
bi <- fread("D://Dropbox/FS/SourceData/BreedingDB__CHECK.csv")
bi$Class <- as.factor(bi$Class)
bi$long <- bi$Lon
bi$lat <- bi$Lat

bat.symbol <- c.img[[1]]
cols <-  c("green2", "dodgerblue2", "magenta1")
col.list <- list()
for (i in 1:nrow(bi)) {
  j <- match(bi$Class[i], levels(bi$Class))
  col.list[[i]] <- cols[[j]]
}
##Sub frames for smooth

bi.tbl <- bi %>%
  group_by(Class, Start) %>%
  summarise(num = n())
bi.ptr <- filter(bi.tbl, Class == 1)
bi.mic <- filter(bi.tbl, Class == 3)
bi.mol <- filter(bi.tbl, Class == 2)

bi.t <- bi %>%
  group_by(Start) %>%
  summarise(num = n())

#### Pannel 1 Map ####
bi.plot <- ggplot()+ bkg +
  geom_polygon(data = fortify(afr.poly),
               aes(long, lat, group = group), 
               colour = "grey20",
               alpha = .25) +
  geom_polygon(data = fortify(rf.poly),
               aes(long, lat, group = group),
               colour = "white", 
               alpha = .25,
               fill = "cornsilk")+
  coord_fixed(xlim = c(-20, 53),ylim = c(-36, 40))

for(i in 1:nrow(bi)){
  bi.plot <- bi.plot +
    add_phylopic(bat.symbol, 1,
                 x = bi$long[i],
                 y = bi$lat[i],
                 ysize = 1.7, 
                 color = col.list[[i]])
}

bi.plot ## looks pretty good. still need to sort the ledgends though

#### Pannel 2 Bar with Smooth ####
bi.bar <- ggplot(data = bi)+ geom_bar(aes(x=Start, fill = Class)) +
  # geom_smooth(data = bi.ptr, aes(x=Start, y = num, color = Class)) #+
  # geom_smooth(data = bi.mic, aes(x=Start, y = num, color = Class)) #+
  # geom_smooth(data = bi.mol, aes(x=Start, y = num, color = Class)) 
  geom_smooth(data = bi.t, aes(x=Start, y = num))
### Need to figure out this thing... total looks alright but broken by class it has 
### huge margins


#### Figure 3 ####
# left column: 3 maps with the annual force of birthing surfaces
# right column: stratafied seasonality smooths (pulses over the year within 5 classes?)
# figures should be based on the *.dbl.imp layers I reckon...

sumGen <- function(model.string){
  ## function for loading rasters and producing an averaged product based on the
  ## model string argument
  f.list <- file.path(clean.dir,"BirthForce",list.files(file.path(clean.dir,"BirthForce"),
                                           pattern = paste0("BR_", model.string)))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again")
  }
  stk <- stack(f.list)
  m.stk <- mean(stk)
}
ptr.stk <- sumGen("ptr.dbl.imp")
mol.stk <- sumGen("mol.dbl.imp")
mic.stk <- sumGen("mic.dbl.imp")

## is this what we want? they are really smoothed thanks to the log transformation...