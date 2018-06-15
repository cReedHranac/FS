#### Having fun with the data
source("R/helperFunctions.R")

#### Packages #####
library(tidyr); library(dplyr); library(raster); library(rgdal); library(ggplot2)
library(gtools); library(data.table)

#### functions ####
better.names <- function(x){
  ### function for impoving names accociated with items retrieved from SpatHandler
  base <- substring(names(x[[1]]), 1, 3)
  i <- 1:12
  j <- c(i[12],i[1:11])
  names(x) <- paste0(base, "_", j, "_", i)
  return(x)
}

spatHandler <- function(model.string){
  ## function for loading rasters produced from spatGLM and producing an averaged product based on the
  ## model string argument
  f.list <- file.path(mod.out.dir, "spatGLM", list.files(file.path(mod.out.dir,"spatGLM"),
                                                         pattern = paste0(model.string,"_")))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again \n *cough* dumbass *cough*")
  }
  ## order and read
  o.list <- mixedsort(f.list)
  stk <- stack(o.list)
  stk <- better.names(stk)
  m.stk <- mean(stk)
  out.l <- list(stk,m.stk)
  return(out.l)
}

rankRisk <- function(x){
  ## function to create risk ranking from spatHandler objects and prep objects for plotting 
  
  ##create naming convention
  base <- substring(names(x[[1]]), 1, 3)[[1]]
  
  ##create dataframes
  x.mean <- data.frame(rasterToPoints(x[[2]], spatial = T))
  colnames(x.mean)[2:3] <- c("long", "lat")
  x.all <- cbind(data.frame(rasterToPoints(x[[1]], spatial = T)),
                 cell = cellFromXY(object = x[[2]],
                                   xy = rasterToPoints(x[[1]], spatial = T)))
  colnames(x.all)[1:2] <- c("long", "lat")
  
  ##mean rank
  out.mean <- x.mean %>% 
    mutate(pct.rank = percent_rank(layer)) %>%
    dplyr::select(-optional)
  
  
  ##All rank
  out.all <- x.all %>%  
    tidyr::gather(key = "window", value = "Rank",
            starts_with(base), factor_key = T) %>%
    mutate(pct.rank = percent_rank(Rank))%>%
    dplyr::select(-optional)
  
  return(list(out.mean, out.all))
}

#### data ####
## read in human spatGLM results
hum.ob <- spatHandler("humNoAn")

hum.rank <- rankRisk(hum.ob)

afr.poly <- readOGR(dsn = path.expand(file.path(data.source, "Africa")),
                    layer = "AfricanCountires")

#### plots ###

ggplot() +
  aes(x= long, y= lat)+
  geom_raster(data = hum.rank[[1]],
              aes(fill = pct.rank))+
  geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "black",
                 fill = NA,
                 size=0.1) +
  coord_equal() +
  theme_bw()
  




Africa.ext <- c(-18, 47, -36, 16)

drc.hd <- readOGR(file.path(data.source, "zone_ste_puc"),
                  "Zone_Ste_Puc")
## proj4 I use
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
drc <- spTransform(drc.hd, CRS(wgs))
rm(drc.hd)
## human outbreak risk w/o animal stocasticity
## animals outbreak risk
ann.ob <- spatHandler("ann")

h.ob.df <- as.data.frame(hum.ob[[1]], row.names = 1:ncell(hum.ob[[1]]))
h.ob.df$cell <- rownames(h.ob.df)



