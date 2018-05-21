#### May DRC EBV outbreak analysis ####

source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

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

p1 <- function(x, y = drc, province = NULL, district = NULL){
  ## function for plotting the DRC health boards with the Risk plots (average)
  ## x <- results from spatHandler
  ## y <- drc
  ## level <- variable to determine the level we're looking at
  ## province <- option for province level selection
  ## distric <- option for province level selection
  
  ## Plot level selection
  if(is.null(province) && is.null(district)){
    y.level <- y
  } 
  if(!is.null(province) && is.null(district)){
    y.level <- y[y@data$PROVINCE == province,]
  }
  if(!is.null(province) && !is.null(district)){
    y.level <- y[y@data$PROVINCE == province & y@data$Nom_ZS_PUC == district,]
  }
  
  ## Manage raster
  x.crop <- mask(crop(x[[2]], extent(y.level)), y.level)
  x.spatial <- rasterToPoints(x.crop, spatial = T)
  y.level$Avg.Risk <- over(y.level, x.spatial, fn = mean)
  y.level$Bikoro <- y.level$Nom_ZS_PUC =="Bikoro" #add bit for our zone of intrest
  
  ## manipulate data structure
  y.level$id <- rownames(y.level@data)
  y.pts <- fortify(y.level, region = "id")
  y.df <- plyr::join(y.pts, y.level@data, by = "id")
  y.df$Risk <- y.df$Avg.Risk$layer
  y.rank <- y.df %>%
    dplyr::select(-Avg.Risk) %>%
    mutate(rank = percent_rank(Risk))
  
  ##plot
  p <- ggplot(y.rank) +
    aes(long,lat,group = group, fill = rank, color = Bikoro) + 
    geom_polygon() +
    coord_fixed() +
    scale_fill_gradient2(low = "yellow", high = "red4",
                         # limits = c(1e-4, 12),
                         na.value = "yellow", 
                         name = "Ebola \nSpillover \nRisk",
                         guide= "colourbar")+
    theme_bw()
  
  return(p)
}

p2 <- function(x, y = drc, province = NULL, district = NULL){
  ## function for plotting the DRC health boards with the Risk plots Monthly
  ## x <- results from spatHandler
  ## y <- drc
  ## level <- variable to determine the level we're looking at
  ## province <- option for province level selection
  ## distric <- option for province level selection
  
  OB.ZONE <- y$Nom_ZS_PUC =="Bikoro" #add bit for our zone of intrest
  y.z <- maptools::spCbind(y, OB.ZONE) 
  
  ## Plot level selection
  if(is.null(province) && is.null(district)){
    y.level <- y.z
  } 
  if(!is.null(province) && is.null(district)){
    y.level <- y.z[y.z@data$PROVINCE == province,]
  }
  if(!is.null(province) && !is.null(district)){
    y.level <- y.z[y.z@data$PROVINCE == province & y.z@data$Nom_ZS_PUC == district,]
  }
  ## Creating names 
  base <- substring(names(x[[1]]), 1, 3)
  i <- 1:12
  j <- c(i[12],i[1:11])
  names(x[[1]]) <- paste0(base, "_", j, "_", i)
  
  ## Manage raster
  x.crop <- mask(crop(x[[1]], extent(y.level)), y.level)
  x.spatial <- rasterToPoints(x.crop, spatial = T)
  Risk.dat <- over(y.level, x.spatial, fn = mean)
  y.dat <- maptools::spCbind(y.level,Risk.dat)
  
  ## manipulate data structure to df
  y.dat$id <- rownames(y.dat@data)
  y.pts <- fortify(y.dat, region = "id")
  y.df <- plyr::join(y.pts, y.dat@data, by = "id")
  
  ## create ranks
  y.rank <- y.df %>%
    dplyr::select(starts_with(substring(names(x.crop),1,3)[[1]])) %>%
    mutate_all(percent_rank)
  colnames(y.rank) <- paste0(base, "_", j, "_", i, "_", "Rank")
  y.full <- cbind(y.df,y.rank)
  
  ## modify for facet plotting
  res.long <- tidyr::gather(data = y.full,
                            key = "window", value = "Risk",
                            ends_with("Rank"), factor_key = T)
  
  ##plot
  p <- ggplot(res.long) +
    aes(long,lat,group = group, fill = Risk, color = OB.ZONE) + 
    geom_polygon() +
    coord_fixed() +
    scale_fill_gradient2(low = "yellow", high = "red4",
                         # limits = c(1e-4, 12),
                         na.value = "yellow", 
                         name = "Ebola \nSpillover \nRisk") + 
    theme_bw() +
    theme( axis.title = element_blank()) +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~ window, ncol = 3)
  
  return(p)
}

p3 <- function(x, region.extent= Africa.ext, afr = afr.poly, drc.hb = drc){
  ###Function for plotting new DRC outbreak against entire africa (average)
  if(is.null(region.extent)){
    x.in <- x[[2]]
  }else{
    x.in <- crop(x[[2]],region.extent)
  }
  ##something stange occurs if crop has be run, appears NA are given data values?
  x.rnk <- calc(x.in, fun = function(x) rank(x, na.last = "keep"))
  x.df <- data.frame(rasterToPoints(x.rnk))
  colnames(x.df) <- c("long","lat","Rank")
  
  
  ##select just the single zone for the drc
  dz1 <- drc.hb[drc.hb$Nom_ZS_PUC == "Bikoro",]
  
  ##Plot
  g.rank <- ggplot(x.rank, aes(x = long, y = lat))+
    # background
    geom_polygon(data = fortify(afr),
                 aes(group = group), 
                 colour = NA,
                 fill = 'black',
                 alpha = .2) +
    
    # fill data values
    geom_raster(aes(fill = Rank), interpolate=TRUE) +
    
    # color scale
    scale_fill_gradient2(low = "yellow", high = "red4",
                         na.value = "yellow", 
                         name = "Mean Rank \nEbola Risk") +
    # country borders
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "black",
                 fill = NA,
                 size=0.1) +
    
    #borders for OB zone
    geom_polygon(data = fortify(dz1),
                 aes(long, lat, group = group), 
                 colour = "green4",
                 fill = NA,
                 size=0.1) +
    
    # limit coordinates
    coord_fixed(xlim = region.extent[1:2], ylim = region.extent[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    # theme
    theme_bw() +
    theme(axis.title = element_blank())
    
  return(g.rank)
}

p4 <- function(x, region.extent= Africa.ext, afr = afr.poly, drc.hb = drc){
  ###Function for plotting new DRC outbreak against entire africa (average)
  base <- substring(names(x[[1]][[1]]),1,3)
  
  if(is.null(region.extent)){
    x.in <- x[[1]]
  }else{
    x.in <- crop(x[[1]],region.extent)
  }
  ##something stange occurs if crop has be run, appears NA are given data values?
  x.df <- data.frame(rasterToPoints(x.in))
  colnames(x.df)[1:2] <- c("long","lat")
  x.rnk <- x.df %>%
    dplyr::select(starts_with(base)) %>%
    mutate_all(rank, na.last = "keep")
  x.full <- cbind(x.df[,1:2], x.rnk)
  
  
  ## modify for facet plotting
  res.long <- tidyr::gather(data = x.full,
                            key = "window", value = "Rank",
                            starts_with(base), factor_key = T)
  
  
  ##select just the single zone for the drc
  dz1 <- drc.hb[drc.hb$Nom_ZS_PUC == "Bikoro",]
  
  ##Plot
  g.rank <- ggplot(res.long, aes(x = long, y = lat))+
    
    # background
    geom_polygon(data = fortify(afr),
                 aes(group = group), 
                 colour = NA,
                 fill = 'black',
                 alpha = .2) +
    
    # fill data values
    geom_raster(aes(fill = Rank), interpolate=TRUE) +
    
    # color scale
    scale_fill_gradient2(low = "yellow", high = "red4",
                         na.value = "yellow", 
                         name = "Mean Rank \nEbola Risk") +
    # country borders
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "black",
                 fill = NA,
                 size=0.1) +
    
    #borders for OB zone
    geom_polygon(data = fortify(dz1),
                 aes(long, lat, group = group), 
                 colour = "green4",
                 fill = NA,
                 size=0.1) +
    
    # limit coordinates
    coord_fixed(xlim = region.extent[1:2], ylim = region.extent[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    # theme
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap(~window, ncol = 3)
  
  return(g.rank)
}

#### data ####
##DRC health districs
afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")

Africa.ext <- c(-18, 47, -36, 16)

drc.hd <- readOGR(file.path(data.source, "zone_ste_puc"),
                  "Zone_Ste_Puc")
## proj4 I use
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
drc <- spTransform(drc.hd, CRS(wgs))
rm(drc.hd)
## human outbreak risk w/o animal stocasticity
hum.ob <- spatHandler("humNoAn")
## animals outbreak risk
ann.ob <- spatHandler("ann")

#### Exploration ####

## we know it's in Equateur provice Bikoro

ob.zone <- drc[drc$PROVINCE == "Equateur" & drc$Nom_ZS_PUC == "Bikoro",]

zone.cells <- as.data.frame(cellFromPolygon(hum.ob[[2]], ob.zone, weights = T))

##Make better names for outbreak


## Convert to df
h.ob.df <- as.data.frame(hum.ob[[1]], row.names = 1:ncell(hum.ob[[1]]))
h.ob.df$cell <- rownames(h.ob.df)


res.long <- h.ob.df %>%
  tidyr::gather(key = "window", value = "Rank",
                starts_with(base), factor_key = T) %>%
  mutate(pct.rank = percent_rank(Rank)) %>%
  # dplyr::filter(cell %in% zone.cells$cell,
  #               window %in% c("hum_3_4", "hum_4_5")) %>%
  dplyr::select(window, pct.rank)

write.csv(res.long, "data/DRC_OB_May2018_Rank.csv", row.names = F)


## Lets plot
p1(x = hum.ob, y = drc)
p2(hum.ob, y = drc)
p3(hum.ob)
p4(hum.ob)

