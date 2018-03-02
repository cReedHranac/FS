###############################################################################
#### Figures for publication  29/01/2018                               ########
###############################################################################
source("R/helperFunctions.R")

library(ggplot2); library(dplyr); library(data.table); library(gridExtra)
#### Figure 1 ####
#3 pannel figure, 1 pannel, map w/ outbreak locations (symbols for an/hum), 
# and 2 outbreak time lines, 1 entire history, one collapsed year
  #### data ####
an.ob <- fread(file.path(clean.dir, "annOB.PPM.csv")) #animal points
an.an <- fread(file.path(data.source, "vIndex", "Animal_Index_Jan2018.csv")) #notes on outbreaks
an.sub <- an.an %>%
  dplyr::select(OUTBREAK_ID, Org.smp, Year.Start)
names(an.sub)[1] <- names(an.ob)[1]
an.full <- an.sub %>%
  inner_join(an.ob,an.sub, by = "Outbreak_ID")

hum.ob <- fread(file.path(clean.dir, "humOB.PPM.csv"))
hum.an <- fread(file.path(data.source, "vIndex","Human_Index_12_2_18.csv"))
hum.sub <- hum.an %>%
  dplyr::select(Outbreak_ID, Year.Start)
hum.full <- hum.sub %>%
  inner_join(hum.ob, hum.sub, by = "Outbreak_ID")


hum.full$Org.smp <- "human"
hum.s <- hum.full %>%
  dplyr::select(which(names(hum.full) %in% names(an.full)))

ob.full <- bind_rows(hum.s, an.full)
colnames(ob.full)[3:4] <- c( "long", "lat")
ob.full$Org.smp <- as.factor(ob.full$Org.smp)

z <- list()
for(i in 1:nrow(ob.full)){
  ifelse(ob.full[i,"Org.smp"] == "human", z[[i]] <- "human",z[[i]] <-  "animal")  
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
  j <- match(ob.full[i,"Org.smp"], levels(ob.full$Org.smp)) 
  img.list[[i]] <- c.img[[j]]
}

## create color vector
cz <- c("darkorange2", "black", "green3", "dodgerblue2", "magenta1")
col.list <- list()
for(i in 1:nrow(ob.full)){
  j <- match(ob.full[i,"Org.smp"], levels(ob.full$Org.smp)) 
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
  coord_fixed(xlim = c(-20, 53),ylim = c(-36, 15))
  
  
for( i in 1:nrow(ob.full)){
  if(ob.full$Org.smp[i] == "human"){
    j <- 3
  } else{ifelse(ob.full$Org.smp[i] == "bat",  j <- 1.9,  j <- 2.5)}
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
  if(ob.full$Org.smp[i] == "human"){
    j <- .3
  } else{ifelse(ob.full$Org.smp[i] == "bat",  j <- .2,  j <- .25)}
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
ob.T$Date <- as.yearmon(paste(ob.T$Year.Start, ob.T$Month.Start), "%Y %m")
ob.a <- ob.T %>%
  dplyr::arrange(Date)


## Plot 

g.time <- ggplot(data = ob.a, aes(x= Date, y = 0))+
  geom_point(data = dplyr::select(ob.a, Date), ## Grey points
             aes(x = Date, y= 0, size = 1.25, alpha = .5),
             color = "grey70", shape=20,
             show.legend = F)+
  geom_point(aes(x = Date, y= 0, color = Org.smp, alpha = .5, size = 1.25),
             show.legend = F)+
  geom_segment(aes(x = 1975, y = 0, xend = 2018, yend = 0),
                arrow = arrow(length =  unit(x = 0.2,units = 'cm'),type = 'closed')) +
  scale_x_yearmon(format = "%Y", n = 10) +
  facet_wrap(~Org.smp, ncol = 1) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank())

g.time
#### Pannel 3 Bar with Smooth ####
### If we can start this at setptember we can likely get the bimodal distribution
### to show up better. Do two different smooths for human/nonhumans and maybe one for total?
ob.hist <- ob.a

org=data.frame(Org.smp=levels(ob.hist$Org.smp),
               Org=c("Chiroptera",
                       "Non-Volant Mammals",
                       "Non-Volant Mammals",
                       "Non-Volant Mammals",
                       "Human"))
ob.hist <- inner_join(ob.hist, org, "Org.smp")

ob.hist$Month <- factor(format(ob.hist$Date, "%b"), 
                        levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

ggplot(data=ob.hist)+ geom_bar(aes(x=Month,fill=Org.smp))+
  facet_wrap(~Org, nrow = 3) +
    bkg + 
  theme(
    axis.ticks.x = element_blank()
  )

#### Figure 2 ####
# Map pannel(s) with bat locations, pannel with collaped year joyplots

  ### Data #### 
bi <- fread("data/afrBatBirthDB.csv")
bi$Class <- factor(bi$Class,
                      labels= c("Pteropodiae", "Molosidae", "Non-Molosid Microbats"))
bi$Start <- factor(bi$Start,
                   labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
bi$long <- bi$Lon
bi$lat <- bi$Lat

bat.symbol <- c.img[[1]]
cols <-  c("green2", "dodgerblue2", "magenta1")
col.list <- list()
for (i in 1:nrow(bi)) {
  j <- match(bi$Class[i], levels(bi$Class))
  col.list[[i]] <- cols[[j]]
}


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
    add_phylopic(bat.symbol, .65,
                 x = bi$long[i],
                 y = bi$lat[i],
                 ysize = 1.7, 
                 color = col.list[[i]])
}

bi.plot ## looks pretty good. still need to sort the ledgends though

#### Pannel 2 Bar with Smooth ####
bi.bar <- ggplot(data = bi)+
  geom_bar(aes(x=Start, fill = Class),
           show.legend = F) +
  facet_wrap(~Class, ncol = 1, scales = "free") +
  bkg +
  theme(
    axis.ticks.x = element_blank())
  
bi.bar
#### Figure 3 BrithForce Model Results####
# left column: 3 maps with the annual force of birthing surfaces
# right column: stratafied seasonality smooths (pulses over the year within 5 classes?)
# figures should be based on the *.dbl.imp layers I reckon...

#### data ####
sumGen <- function(model.string){
  ## function for loading rasters and producing an averaged product based on the
  ## model string argument
  f.list <- file.path(clean.dir,"BirthForce",list.files(file.path(clean.dir,"BirthForce"),
                                           pattern = paste0("BR_", model.string)))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again \n *cough* dumbass *cough*")
  }
  ## stupid hack to order list since names are too complex for mixed sort
  o.list <- f.list[c(1,5,6,7,8,9,10,11,12,2,3,4)]
  stk <- stack(o.list)
  m.stk <- mean(stk)
  out.l <- list(stk,m.stk)
  return(out.l)
}
ptr.sum <- sumGen("ptr.dbl.imp")
mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")

#### Pannel 1 (Left) ####

BFgplot <- function(x, source.path = data.source){
  #### Set Up ####
  ## accessory layers
  afr.poly <- readOGR(dsn = file.path(source.path, "Africa"),
                      layer = "AfricanCountires")
  rf.poly <- rasterToPolygons(raster(file.path(source.path, "cropMask.tif")),
                              fun = function(x){x==1}, dissolve = T)
  
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Number")
  
  g.plot <- ggplot(sum.df) +
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 alpha = .25) +
    aes(x=long, y=lat) +
    scale_fill_gradient(low = "yellow", high = "red4",
                        limits = c(0,max(sum.df$Number)))+
    geom_raster(aes(fill = Number), interpolate = T)+
    
    #add area modeled
    geom_polygon(data = fortify(rf.poly),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 fill = NA,
                 alpha = .2) +
    coord_fixed(xlim = c(-20, 53),ylim = c(-36, 15))
  
}

ptr.BF <- BFgplot(x = ptr.sum)
mol.BF <- BFgplot(x = mol.sum)
mic.BF <- BFgplot(x = mic.sum)

#### Pannel 2 (Right) ####
# install_github("cran/ggridges")
library(ggridges); library(readr)
sub.ext <- c(-20, 53, -36, 15) #extent subset like that of the other map figures
BFridge <- function(x, n.bin, crop.extent = sub.ext){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.crop <- crop(x[[1]], crop.extent)
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  colnames(x.cv)[3:ncol(x.cv)] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  x.df <- x.cv[complete.cases(x.cv),]
  x.df$Jan2 <- x.df$Jan
  bf.df <- x.df %>%
    gather("month","BF",3:ncol(x.df)) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(bf.mean = mean(BF)) 

   bf.df$month <-  factor(bf.df$month,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec", "Jan2"))
  
  
  bf.ridge <- ggplot(data= bf.df, 
                     aes(x= month,y= strata,height = bf.mean, group = strata, fill = bf.mean))+
    geom_density_ridges_gradient(stat = "identity", scale = 3, alpha = .5) +
    scale_fill_gradient(low = "yellow", high = "red4",
                        limits = c(0,max(bf.df$bf.mean))) +
    scale_x_discrete(label = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
    
    
    
    bkg <- theme(
      panel.background = element_rect(fill = "lightblue",
                                      colour = "lightblue",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = "white"),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text( hjust = .5))
  return(bf.ridge + bkg)
}

r.ptr <- BFridge(x = ptr.sum, n.bin = 40, crop.extent = sub.ext)
r.mic <- BFridge(mic.sum, 40)
r.mol <- BFridge(mol.sum, 40)


#### Figure 4 SpatGLM Model results ####
#Same as birthForce model but obviously, human spatGLM results

#### data ####
spatHandler <- function(model.string){
  ## function for loading rasters produced from spatGLM and producing an averaged product based on the
  ## model string argument
  f.list <- file.path(mod.out.dir, "spatGLM", list.files(file.path(mod.out.dir,"spatGLM"),
                                                        pattern = paste0(model.string)))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again \n *cough* dumbass *cough*")
  }
  ## stupid hack to order list since names are too complex for mixed sort
  o.list <- mixedsort(f.list)
  stk <- stack(o.list)
  m.stk <- mean(stk)
  out.l <- list(stk,m.stk)
  return(out.l)
}

#### Pannel 1 Average ####
bkg <- theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  plot.title = element_text(hjust = 0.5))
mylog = scales::trans_new('mylog',
                          transform=function(x) { log(x+1e-5) },
                          inverse=function(x) { exp(x)-1e-5})

ERgplot <- function(x, source.path = data.source){
  #### Set Up ####
  ## accessory layers
  afr.poly <- readOGR(dsn = file.path(source.path, "Africa"),
                      layer = "AfricanCountires")
  rf.poly <- rasterToPolygons(raster(file.path(source.path, "cropMask.tif")),
                              fun = function(x){x==1}, dissolve = T)
  
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Risk")
  sum.df$Risk <- round(sum.df$Risk, 3)
  g.plot <- ggplot(sum.df) +
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 alpha = .25) +
    aes(x=long, y=lat) +
    geom_raster(aes(fill = Risk), interpolate = T)+
    scale_fill_gradientn(name = "Risk", trans = mylog,
                         colors = terrain.colors(10)) +
                        # low = "yellow", high = "red4",
                        # breaks = c(1e-4,.01, .05, .07, .1))+
    #add area modeled
    geom_polygon(data = fortify(rf.poly),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 fill = NA,
                 alpha = .2) +
    coord_fixed(xlim = c(-20, 53),ylim = c(-36, 15))
  out <- g.plot + bkg
}

hum.mean <- spatHandler("hum")
risk.plot <- ERgplot(hum.mean)
risk.plot
#### Pannel 2 (Right) ####
library(ggridges)
sub.ext <- c(-20, 53, -36, 15) #extent subset like that of the other map figures
ERridge <- function(x, n.bin, crop.extent = sub.ext){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.crop <- crop(x[[1]], crop.extent)
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  colnames(x.cv)[3:ncol(x.cv)] <- c("January","February","March",
                                    "April","May","June","July","August","September",
                                    "October","November","December")
  x.df <- x.cv[complete.cases(x.cv),]
  ER.df <- x.df %>%
    gather("month","ER",3:14) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(ER.mean = mean(ER)) 
  
  ER.df$month <-  factor(ER.df$month,levels=c("January","February","March",
                                              "April","May","June","July","August","September",
                                              "October","November","December"))
  
  
  ER.ridge <- ggplot(data= ER.df, 
                     aes(x= month,y= strata,height = ER.mean, group = strata, fill = ER.mean))+
    geom_density_ridges_gradient(stat = "identity", scale = 2.5) +
    # scale_fill_gradient(name = "Risk", trans = "log10",
    #                     low = "yellow", high = "red4",
    #                     breaks = c(0,.001,.01))
    scale_fill_gradient(low = "yellow", high = "red4",
                        limits = c(0,max(ER.df$ER.mean)))

  
  
  bkg <- theme(
    panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1))
  return(ER.ridge + bkg)
}
k <- ERridge(hum.mean, n.bin = 40, sub.ext)
k
