##########################
#### Figure 4         ####
##########################
source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges)
#### functions ####
spatHandler <- function(model.string){
  ## function for loading rasters produced from spatGLM and producing an averaged product based on the
  ## model string argument
  f.list <- file.path(mod.out.dir, "spatGLM", list.files(file.path(mod.out.dir,"spatGLM"),
                                                         pattern = paste0(model.string,"_")))
  
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

scaleFUN <- function(x)sprintf("%.4f", round(x, digits = 4)) 

ERgplot <- function(x, source.path = data.source, afr = afr.poly, rf = rf.poly){
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Risk")
  
  g.plot <- ggplot(sum.df) +
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 alpha = .2) +
    aes(x=long, y=lat) +
    
    #Color
    scale_fill_gradient2(trans = scales::log_trans(),
                        low = "yellow", mid = "red4",
                        #limits = c(1e-20, 2),
                        limits = c(1e-4, 12),
                        na.value = "yellow", 
                        name = "Mean \nEbola \nRisk") +
    
    #fill data values
    geom_raster(aes(fill = Risk), interpolate = T)+
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 fill = NA,
                 alpha = .2) +
    
    #add area modeled
    geom_polygon(data = fortify(rf),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    #Extras
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15)) +
    scale_y_continuous(expand = c(0,0) )#, lables = scaleFUN) +
  theme_bw()+
    theme(axis.title = element_blank())
  
  coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15))
  out <- g.plot + theme_bw()
}

ERridge <- function(x, n.bin, crop.extent = sub.ext){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.crop <- raster::crop(x[[1]],y = raster::extent(crop.extent))  
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  colnames(x.cv)[3:ncol(x.cv)] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  x.df <- x.cv[complete.cases(x.cv),]
  x.df$Jan2 <- x.df$Jan
  ER.df <- x.df %>%
    tidyr::gather("month","ER",3:ncol(x.df)) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(ER.mean = mean(ER))
  
  ER.df$month <-   factor(ER.df$month,levels=c("Jan","Feb","Mar",
                                               "Apr","May","Jun",
                                               "Jul","Aug","Sep",
                                               "Oct","Nov","Dec", "Jan2"))
  
  ER.ridge <- ggplot(data= ER.df,
                     aes(x= month,y= strata,height = ER.mean, group = strata, fill = ER.mean))+
    geom_density_ridges_gradient(stat = "identity", scale = 5) +
    scale_fill_gradient(low= "yellow", high = "red4",
                        limits = c(0,max(ER.df$ER.mean)),
                        name = "Mean \nEbola \nRisk") +
    scale_x_discrete(label = c("Jan","Feb","Mar",
                               "Apr","May","Jun",
                               "Jul","Aug","Sep",
                               "Oct","Nov","Dec", "Jan"),
                     expand = c(0,0))+
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
  return(ER.ridge)
}
#### Pannel 1 Average ####

afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")
rf.poly <- rasterToPolygons(raster(file.path(data.source, "cropMask.tif")),
                            fun = function(x){x==1}, dissolve = T)

hum.mean <- spatHandler("hum") #This for averages 
risk.plot <- ERgplot(hum.mean)

### Subregion plots ###
  #Central
central.africa <- c(8,35,-5,6)
CentralZone <- risk.plot + coord_fixed(xlim = central.africa[1:2],ylim = central.africa[3:4]) +
  theme(legend.position="none") + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

  #Western
west.africa <- c(-14,5, 4, 12)
westernZone <- risk.plot + coord_fixed(xlim = west.africa[1:2],ylim = west.africa[3:4]) +
  theme(legend.position="none") + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"))


#### Alternative human models for ridges ####

hum.NoAn <- spatHandler("humNoAn") # Go with this one for ridges

Afr.ext <- c(-18, 49, -36, 15)
afr.ridge <- ERridge(hum.NoAn, n.bin = 50, crop.extent = Afr.ext )
central.ridge <- ERridge(hum.NoAn, n.bin = 30, central.africa)
western.ridge <- ERridge(hum.NoAn, n.bin = 30, west.africa)

