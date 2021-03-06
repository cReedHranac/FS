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

spatHandler <- function(model.string, mod.dir){
  ## function for loading rasters produced from spatGLM and producing an averaged product based on the
  ## model string argument
  base <- substring(model.string, 1, 3)
  
  f.list <- mixedsort(list.files(file.path(mod.out.dir,mod.dir),
                                 pattern = paste0(model.string,"_"), 
                                 full.names = T))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again \n *cough* dumbass *cough*")
  }
  ## order and read
  stk <- better.names(stack(f.list))
  m.stk <- mean(stk)
  out.l <- list(stk,m.stk)
  return(out.l)
}

p1 <- function(x, y = drc, outbreaks = NULL, province = NULL, district = NULL){
  ## function for plotting the DRC health boards with the Risk plots (average)
  ## x <- results from spatHandler
  ## y <- drc
  ## outbreaks <- list of zones to highlght
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
  y.level$Outbreaks <- y.level$Nom_ZS_PUC %in% c(outbreaks) #add bit for our zone of intrest
  
  ## manipulate data structure
  y.level$id <- rownames(y.level@data)
  y.pts <- fortify(y.level, region = "id")
  y.df <- plyr::join(y.pts, y.level@data, by = "id")
  y.df$Risk <- y.df$Avg.Risk$layer
  y.rank <- y.df %>%
    dplyr::select(-Avg.Risk) %>%
    mutate(rank = (percent_rank(Risk)*100))
  y.rank$Outbreaks[which(y.rank$Outbreaks != F)] <- as.character(y.rank$Nom_ZS_PUC[which(y.rank$Outbreaks != F)])
  y.rank$Outbreaks[which(y.rank$Outbreaks == F)] <- "None"
  y.rank$Outbreaks <- as.factor(y.rank$Outbreaks)
  
  ##create unify for border
  y.simple <- rgeos::gSimplify(y, tol = .0000001)
  y.buff <- rgeos::gBuffer(y.simple, width=0)
  y.uni <- rgeos::gUnaryUnion(y, id = "Province")
  
  ##plot
  p <- ggplot(y.rank) +
    ##poly fill
    aes(long,lat,group = group, fill = rank) + 
    ##poly lines
    geom_polygon(aes(color = Outbreaks), size = 1) +
    guides(color = "none")+ #kill lables
    # geom_line(data = data.frame(x = c(rep(min(y.rank$long),100)),
    #                             y = c(seq(0,1, length.out = 100))),
    #                             aes(x=x, y=y),
    #                             size = 2)+
    coord_fixed() +
    ##color for gradient
    scale_fill_gradient2(low = "yellow", high = "red4",
                         limits = c(0,100),
                         na.value = "white", 
                         name = "EVD \nSpillover \nPercent\n Rank",
                         guide= "colourbar")+
    ##color for poly lines (last is transparent)
    scale_color_manual(values=c("#E69F00", "#56B4E9", "#ffffff00")) +
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = c(0.11,0.76),
          legend.margin = margin(),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.4, "cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=9),
          plot.margin=unit(rep(0,4), "cm"))

  
  return(p)
}

p1raw <- function(x, y = drc, outbreaks = NULL, province = NULL, district = NULL){
  ## function for plotting the DRC health boards with the Risk plots (average)
  ## x <- results from spatHandler
  ## y <- drc
  ## outbreaks <- list of zones to highlght
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
  y.level$Outbreaks <- y.level$Nom_ZS_PUC %in% c(outbreaks) #add bit for our zone of intrest
  
  ## manipulate data structure
  y.level$id <- rownames(y.level@data)
  y.pts <- fortify(y.level, region = "id")
  y.df <- plyr::join(y.pts, y.level@data, by = "id")
  y.df$Risk <- y.df$Avg.Risk$layer
  y.rank <- y.df %>%
    dplyr::select(-Avg.Risk) %>%
    mutate(relRank = (Risk/backgroundGeo))
  y.rank$Outbreaks[which(y.rank$Outbreaks != F)] <- as.character(y.rank$Nom_ZS_PUC[which(y.rank$Outbreaks != F)])
  y.rank$Outbreaks[which(y.rank$Outbreaks == F)] <- "None"
  y.rank$Outbreaks <- as.factor(y.rank$Outbreaks)
  
  ##create unify for border
  y.simple <- rgeos::gSimplify(y, tol = .0000001)
  y.buff <- rgeos::gBuffer(y.simple, width=0)
  y.uni <- rgeos::gUnaryUnion(y, id = "Province")
  
  ##plot
  p <- ggplot(y.rank) +
    ##poly fill
    aes(long,lat,group = group, fill = relRank) + 
    ##poly lines
    geom_polygon(aes(color = Outbreaks), size = 1) +
    guides(color = "none")+ #kill lables
    # geom_line(data = data.frame(x = c(rep(min(y.rank$long),100)),
    #                             y = c(seq(0,1, length.out = 100))),
    #                             aes(x=x, y=y),
    #                             size = 2)+
    coord_fixed() +
    ##color for gradient
    scale_fill_gradient2(low = "yellow", high = "red4",
                         trans = "log10",
                         na.value = "white", 
                         name = "EVD \nSpillover \nRelative\nRisk",
                         guide= "colourbar")+
    ##color for poly lines (last is transparent)
    scale_color_manual(values=c("#E69F00", "#56B4E9", "#ffffff00")) +
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = c(0.11,0.76),
          legend.margin = margin(),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.4, "cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=9),
          plot.margin=unit(rep(0,4), "cm"))
  
  
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
## human outbreak risk w/o animal stocasticity ## running with Prob modification
hum.ob <- spatHandler("humNoAn", "SpGLMRes_Prob")
hum.geo <- calc(hum.ob[[1]], function(x){y <- log(x); return(mean(y))})

backgroundMean <- cellStats(hum.ob[[2]], mean)
backgroundGeo <- exp(cellStats(hum.geo, mean))
## animals outbreak risk
ann.ob <- spatHandler("ann","SpGLMRes_Nov")

#### Exploration ####

## outbreak 1
## we know it's in Equateur provice Bikoro
ob1.zone <- drc[drc$PROVINCE == "Equateur" & drc$Nom_ZS_PUC == "Bikoro",]
zone1.cells <- as.data.frame(cellFromPolygon(hum.ob[[2]], ob1.zone, weights = T))


## outbreak 2
## North Kivu province, zone mangina
ob2.zone <- drc[drc$PROVINCE =="Nord Kivu" & drc$Nom_ZS_PUC == "Beni",]
zone2.cells <- as.data.frame(cellFromPolygon(hum.ob[[2]], ob2.zone, weights = T))

##Make better names for outbreak
## Convert to df
h.ob.df <- as.data.frame(hum.ob[[1]], row.names = 1:ncell(hum.ob[[1]]))
h.ob.df$cell <- rownames(h.ob.df)

base <- "hum"

res.ob1 <- h.ob.df %>%
  tidyr::gather(key = "window", value = "Rank",
                starts_with(base), factor_key = T) %>%
  mutate(pct.rank = percent_rank(Rank)) %>%
  dplyr::filter(cell %in% zone1.cells$cell,
                window %in% c("hum_12_1")) %>%
  dplyr::select(cell, window, pct.rank, Rank)

res.ob1.apr <- h.ob.df %>%
  tidyr::gather(key = "window", value = "Rank",
                starts_with(base), factor_key = T) %>%
  mutate(pct.rank = percent_rank(Rank)) %>%
  dplyr::filter(cell %in% zone1.cells$cell,
                window %in% c("hum_3_4")) %>%
  dplyr::select(cell, window, pct.rank, Rank)


write.csv(res.ob1, "data/DRC_OB_Dec2017_Rank.csv", row.names = F)
write.csv(res.ob1, "data/DRC_OB_Apr2018_Rank.csv", row.names = F)


res.ob2 <- h.ob.df %>%
  tidyr::gather(key = "window", value = "Rank",
                starts_with(base), factor_key = T) %>%
  mutate(pct.rank = percent_rank(Rank)) %>%
  dplyr::filter(cell %in% zone2.cells$cell,
                window %in% c("hum_6_7")) %>%
  dplyr::select(cell, window, pct.rank, Rank)

write.csv(res.ob2, "data/DRC_OB_July2018_Rank.csv", row.names = F)

## Lets plot
ob.rank <- p1(x = hum.ob, y = drc, outbreaks = c("Bikoro", "Beni"))
ob.raw <- p1raw(x = hum.ob, y = drc, outbreaks = c("Bikoro", "Beni"))
# p2(hum.ob, y = drc)
# p3(hum.ob)
# p4(hum.ob)

#### Creating the violin plots ####
## build the data frame
ob1 <- cbind(res.ob1, Outbreak = "Bikoro_Dec", id = paste0(res.ob1$cell,res.ob1$window))
ob1.alt <- cbind(res.ob1.apr, Outbreak = "Bikoro", id = paste0(res.ob1.apr$cell, res.ob1.apr$window))
ob2 <- cbind(res.ob2, Outbreak = "Beni", id = paste0(res.ob2$cell,res.ob2$window))

#mean RR human
hum.mean <- calc(hum.ob[[1]], mean)
hum.mean.df <- as.data.frame(hum.mean, row.names = 1:ncell(hum.mean))
hum.mean.df$cell <- rownames(hum.mean.df)
hum.mean.df$mRR <- percent_rank(hum.mean.df$layer)

## extract cells from each of the outbreaks and ammend to outbreak dataframes
ob1$mRisk <- hum.mean.df$mRR[which(hum.mean.df$cell %in% ob1$cell)]
ob2$mRisk <- hum.mean.df$mRR[which(hum.mean.df$cell %in% ob2$cell)]
ob1.alt$mRisk <- hum.mean.df$mRR[which(hum.mean.df$cell %in% ob1.alt$cell)]
## bind
ob.df <- rbind(ob1.alt,ob2)
ob.df$pct.rank <- ob.df$pct.rank*100


p.vol <- ggplot(ob.df, aes(x = Outbreak, y = pct.rank, color = Outbreak)) + 
  geom_boxplot(size = 1, width = .5) + 
  geom_jitter(shape = 16, position = position_jitter(.1), color = "black")+
  scale_y_continuous(limits = c(0,100)) + 
  labs(x = "Outbreak", 
       y = "Precent Rank")+
  scale_color_manual(values=c("#56B4E9","#E69F00")) + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank())
  

f5.insert <- ob.rank + 
  annotation_custom(grob = ggplotGrob(p.vol),
                    xmin = 11,
                    xmax = 20,
                    ymin = -14.5,
                    ymax = -8.25)


#### Monthly features
#modify names to create continious
h.modify <- h.ob.df
colnames(h.modify) <- c(1:12, "cell")

res.month <- h.modify %>%
   tidyr::gather(key = "window", value = "Rank",
                1:12, factor_key = T, convert = T) %>%
  mutate(pct.rank = percent_rank(Rank)*100) %>%
  dplyr::filter(cell %in% c(zone2.cells$cell, zone1.cells$cell))%>%
  mutate(Outbreak = ifelse(cell %in% zone1.cells$cell, "Bikoro", "Beni")) %>%
  group_by(window, Outbreak) %>%
  mutate(rnk.med = median(pct.rank),
         rnk.low = min(pct.rank),
         rnk.high = max(pct.rank))


(p.rib <- ggplot(data = res.month, aes(x = window, y = pct.rank, colour = Outbreak)) + 
#  geom_point() +
  geom_line(aes(group = cell), alpha = .3) +
  geom_line(aes(y = rnk.med), size=1) +
  scale_x_continuous(breaks = 1:12, labels=substring(month.abb, 1, 1),
                     expand = c(0,0)) +
  ylim(c(0,100))+
  scale_color_manual(values=c("#E69F00","#56B4E9")) +
  guides(color='none') +
    labs(x = "Month", 
       y = "Precent Rank")+
  
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
#        axis.text.x = element_blank(),
        plot.margin = unit(c(0,10,0,0), "points"))
)

## creating the raw data figure
res.raw <-  h.modify %>%
  tidyr::gather(key = "window", value = "Risk",
                1:12, factor_key = T, convert = T) %>%
  mutate(rel.Risk = Risk/backgroundGeo)%>%
  dplyr::filter(cell %in% c(zone2.cells$cell, zone1.cells$cell))%>%
  mutate(Outbreak = ifelse(cell %in% zone1.cells$cell, "Bikoro", "Beni")) %>%
  group_by(window, Outbreak) %>%
  mutate(raw.med = median(rel.Risk),
         raw.low = min(rel.Risk),
         raw.high = max(rel.Risk))

res.raw.ob <-  h.modify %>%
  tidyr::gather(key = "window", value = "Risk",
                1:12, factor_key = T, convert = T) %>%
  mutate(rel.Risk = Risk/backgroundGeo)%>%
  dplyr::filter(cell %in% c(zone2.cells$cell, zone1.cells$cell))%>%
  mutate(Outbreak = ifelse(cell %in% zone1.cells$cell, "Bikoro", "Beni")) %>%
  group_by(Outbreak) %>%
  mutate(raw.med = median(rel.Risk),
         raw.low = min(rel.Risk),
         raw.high = max(rel.Risk))


(p.raw <- ggplot(data = res.raw, aes(x = window, y = rel.Risk, color = Outbreak)) + 
    #  geom_point() +
  geom_line(aes(group = cell), alpha = .3) +
  geom_line(aes(y = raw.med), size=1) +
  # geom_ribbon(aes(ymin = raw.low, ymax = raw.high), alpha = .1, linetype = 2)+
  scale_x_continuous(breaks = 1:12, labels=substring(month.abb, 1, 1),
                     expand = c(0,0)) +
  scale_color_manual(values=c("#E69F00","#56B4E9")) + 
  scale_y_log10()+
  guides(color='none') +
  labs(x = "Month", 
       y = "Relative Risk")+
  theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,10,0,0), "points"))
)
  
  
## Putting the two together
(f5.full <- grid.arrange(ob.rank,p.rib, ob.raw, p.raw,
                        layout_matrix = rbind(c(1,2),
                                              c(3,4)),
                        heights = c(.5,.5)))
plot(f5.full)

ggsave(filename = "figures/Fig5_Prob.pdf",
       plot = f5.full,
       device = cairo_pdf,
       height = 7, 
       width = 7.5,
       units = "in",
       dpi = 300)


ob.df %>% 
  group_by(Outbreak) %>%
  summarise(min = round(min(pct.rank),4), max = max(pct.rank), median = median(pct.rank))

