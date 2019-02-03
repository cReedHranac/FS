#### Functional forms for outbreak investigation and figure 5 generation ####

source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

#### functions ####
better.names <- function(x, model.name){
  ### function for impoving names accociated with items retrieved from SpatHandler
  base <- sapply(strsplit(model.name, "_"), tail, 1)
  i <- 1:12
  j <- c(i[12],i[1:11])
  names(x) <- paste0(base, "_", j, "_", i)
  return(x)
}

spatHandler <- function(model.name, mod.dir){
  ## function for loading rasters produced from spatGLM and producing an averaged product based on the
  ## model string argument
  base <- sapply(strsplit(model.name, "_"), tail, 1)
  
  f.list <- mixedsort(list.files(file.path(mod.dir,model.name),
                                 pattern = "noAn_", 
                                 full.names = T))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again \n *cough* dumbass *cough*")
  }
  ## order and read
  stk <- better.names(stack(f.list), model.name = model.name)
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
  
  geo.mean <- calc(x[[1]],
                   function(x){y <- log(x); return(mean(y))})
  names(geo.mean) <- "geo.mean"
  backgroundGeo <- exp(cellStats(geo.mean, mean))
  
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

fig5.fun <- function(model.name, mod.dir, drc.poly = drc,
                     write.out = T, out.df, out.fig){
  ### Function for generating dataframes and figures associated with the 
  ### original outbreak explorer script (Ob_explore.R)
  ### Fixed for the outbreaks in Birko and Beni
  ## Arguments:
  # model.string <- model string pass to spatHandler
  # mod.dir <- directory to where these models are held
  # drc.poly <- polygonDF including the DHB for the DRC
  # write <- logical flag for if the summary data files should be written out
  # out.df <- where to write dataframe summary to 
  # out.fig <- where the figure should be written out to
  
  ## read in data set
  base <- sapply(strsplit(model.name, "_"), tail, 1)
  
  pred.dat <- spatHandler(model.name = model.name,
                          mod.dir = mod.dir)
  geo.mean <- calc(pred.dat[[1]],
                   function(x){y <- log(x); return(mean(y))})
  names(geo.mean) <- "geo.mean"
  backgroundMean <- cellStats(pred.dat[[2]], mean)
  backgroundGeo <- exp(cellStats(geo.mean, mean))
  
  ## set up outbreaks
  ob.zone1 <- drc.poly[drc.poly$PROVINCE == "Equateur" &
                         drc.poly$Nom_ZS_PUC == "Bikoro",]
  ob.zone2 <- drc.poly[drc.poly$PROVINCE =="Nord Kivu" &
                         drc.poly$Nom_ZS_PUC == "Beni",]
  
  ## define cells
  z1.cells <- as.data.frame(cellFromPolygon(pred.dat[[2]],
                                            ob.zone1, 
                                            weights = T))
  z2.cells <- as.data.frame(cellFromPolygon(pred.dat[[2]],
                                            ob.zone2, 
                                            weights = T))
  
  ## data frame of predicted values
  pred.stk <- stack(pred.dat[[1]],geo.mean)
  pd.df <- as.data.frame(pred.stk,
                         row.names = 1:ncell(pred.dat[[1]]))
  pd.df$cell <- rownames(pd.df)
  names(pd.df)[1:12] <- 1:12
  
  res.month <- pd.df %>%
    tidyr::gather(key = "window", value = "Rank",
                  1:12, factor_key = T, convert = T) %>%
    mutate(pct.rank = percent_rank(Rank)*100) %>%
    dplyr::filter(cell %in% c(z2.cells$cell, z1.cells$cell))%>%
    mutate(Outbreak = ifelse(cell %in% z1.cells$cell, "Bikoro", "Beni")) %>%
    group_by(window, Outbreak) %>%
    mutate(rnk.med = median(pct.rank),
           rnk.low = min(pct.rank),
           rnk.high = max(pct.rank))
  
  res.raw <-  pd.df %>%
    tidyr::gather(key = "window", value = "Risk",
                  1:12, factor_key = T, convert = T) %>%
    mutate(rel.Risk = Risk/backgroundGeo)%>%
    dplyr::filter(cell %in% c(z2.cells$cell, z1.cells$cell))%>%
    mutate(Outbreak = ifelse(cell %in% z1.cells$cell, "Bikoro", "Beni")) %>%
    group_by(window, Outbreak) %>%
    mutate(rel.med = median(rel.Risk),
           rel.low = min(rel.Risk),
           rel.high = max(rel.Risk))
  res.df <- bind_cols(res.month, res.raw[,c("Risk", "rel.Risk",
                                            "rel.low", "rel.med", "rel.high")])
  
  ####Plots####
  (ob.rank <- p1(pred.dat,
                 y = drc.poly,
                 outbreaks = c("Bikoro", "Beni")))
  (ob.raw <- p1raw(pred.dat,
                   y = drc.poly,
                   outbreaks = c("Bikoro", "Beni")))
  (p.rib <- ggplot(data = res.month,
                   aes(x = window,
                       y = pct.rank,
                       colour = Outbreak)) + 
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
            plot.margin = unit(c(0,10,0,0), "points"))
  )
  (p.raw <- ggplot(data = res.raw, aes(x = window, y = rel.Risk, color = Outbreak)) + 
      geom_line(aes(group = cell), alpha = .3) +
      geom_line(aes(y = rel.med), size=1) +
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
  
  #### Write out ####
  (f5.full <- grid.arrange(ob.rank,p.rib, ob.raw, p.raw,
                           layout_matrix = rbind(c(1,2),
                                                 c(3,4)),
                           heights = c(.5,.5)))
  
  if(write.out == T){
    write.csv(res.df,
              file.pat(out.df,
                       paste0("obDF",model.name,".csv")),
              row.names = F)
    ggsave(filename = file.path(out.fig,
                                paste0("Fig5_",model.name,".pdf")),
           plot = f5.full,
           device = cairo_pdf,
           height = 7, 
           width = 7.5,
           units = "in",
           dpi = 300)
  }
  
  out <- list(plot = f5.full,
              p1 = ob.rank,
              p2 = ob.raw,
              p3 = p.rib,
              p4 = p.raw,
              df = res.df)
  return(out)
  
}

altModBoxes.Month <- function(x, df = ob.masterframe,
                              write.out = F, out.fig){
  ### Function for looking at outbreak prediction preformace across alternative models
  ## x <- name of colum from df to use
  ## df <- df to use (default to ob.masterframe)
  ## write.out <- logical, should write
  ## out.fig <- where the figure should be written to
  
  x.enquo <- enquo(x)
  x.quo <- quo(x)
  
  ob.place.time <- df %>%
    filter(Outbreak == "Beni" & window == 6)
  ob.pt <- df %>% 
    filter(Outbreak == "Bikoro" & window == 3) %>%
    bind_rows(ob.place.time)
  
  p.vol <- ggplot(ob.pt,
                  aes(x = Outbreak,
                      y = !!x.enquo,
                      color = Outbreak,
                      fill = Mod.Name)) + 
    geom_boxplot(size = 1, width = .5) + 
    # geom_jitter(shape = 16, 
    #             position = position_jitter(.1),
    #             color = "black")+
#Need to figure out why this doesnt work
    # scale_y_continuous(limits = c(floor(!!x.enquo),
    #                               ceiling(!!x.quo))) + 
    labs(x = "Outbreak", 
         y = quo_name(x.enquo))+
    scale_color_manual(values=c("#56B4E9","#E69F00")) + 
    theme_bw()+
    theme(#legend.position = "none",
      axis.title.x = element_blank())
  
  if(write.out == T){
    a <- file.path(out.fig, paste0("altModBoxesMonth_",quo_name(x.enquo),".pdf"))
    ggsave(a,
           plot = p.vol,
           device = cairo_pdf,
           height = 7,
           width = 7.5,
           units = "in",
           dpi = 300)
  }
  return(p.vol)
  
}
altModBoxes.Total <- function(x, df = ob.masterframe,
                              write.out = F, out.fig){
  ### Function for looking at outbreak prediction preformace across alternative models
  ## x <- name of colum from df to use
  ## df <- df to use (default to ob.masterframe)
  ## write.out <- logical, should write
  x.enquo <- enquo(x)
  x.quo <- quo(x)
  
  p.vol <- ggplot(df,
                  aes(x = Outbreak,
                      y = !!x.enquo,
                      color = Outbreak,
                      fill = Mod.Name)) + 
    geom_boxplot(size = 1, width = .5) + 
    # geom_jitter(shape = 16, 
    #             position = position_jitter(.1),
    #             color = "black")+
    #Need to figure out why this doesnt work
    # scale_y_continuous(limits = c(floor(!!x.enquo),
    #                               ceiling(!!x.quo))) + 
    labs(x = "Outbreak", 
         y = quo_name(x.enquo))+
    scale_color_manual(values=c("#56B4E9","#E69F00")) + 
    theme_bw()+
    theme(#legend.position = "none",
      axis.title.x = element_blank())
  
  if(write.out == T){
    a <- file.path(out.fig, paste0("altModBoxesTotal_",quo_name(x.enquo),".pdf"))
    ggsave(a,
           plot = p.vol,
           device = cairo_pdf,
           height = 7,
           width = 7.5,
           units = "in",
           dpi = 300)
  }
  return(p.vol)
  
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

#### create the figure 5 images ####
human.model.names <- list.files(mod.out.nov,pattern = "h_")

human.fig5s <- lapply(human.model.names, fig5.fun,
                      mod.dir= mod.out.nov,
                      write.out = T,
                      out.df = dOut.1,
                      out.fig = fig.hum1)


## read in the dataframes in and add Model name column
## get names 
dfs <- lapply(list.files(hum.1, pattern = "ob", full.names = T),read.csv)


for(i in 1:length(human.model.names)){
  dfs[[i]]$Mod.Name <- sapply(strsplit(human.model.names[[i]], "_"), tail, 1)
}

ob.masterframe <- do.call(rbind, dfs)
names(ob.masterframe)

t1 <- altModBoxes.Month(pct.rank,
                        write.out = T,
                        out.fig = fig.hum1)
t2 <- altModBoxes.Month(rel.Risk,
                        write.out = T,
                        out.fig = fig.hum1)

t3 <- altModBoxes.Total(pct.rank,
                        write.out = T,
                        out.fig = fig.hum1)
t4 <- altModBoxes.Total(rel.Risk,
                        write.out = T,
                        out.fig = fig.hum1)



### Which model has the hightest median pct.ranks and rel.risk
head(ob.masterframe)

zoop <- ob.masterframe %>%
  filter(Outbreak == "Beni" & window == 6)
zoop2 <- ob.masterframe %>%
  filter(Outbreak == "Bikoro" & window == 3) %>%
  bind_rows(zoop) %>%
  group_by(Mod.Name, Outbreak) %>%
  summarise(m.rank = median(pct.rank),
            m.relRisk = median(rel.Risk)) 
