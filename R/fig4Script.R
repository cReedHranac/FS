##########################
#### Figure 4         ####
##########################

# feel free to specify the path to FS (including the FS folder) in base.path.
# if not set, it'll use what Reed uses (D:\Dropbox\FS or ~/Dropbox/FS)
#base.path <- "~/data/Dropbox/FS"

source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

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

ERridge <- function(x, n.bin, scale = 5, crop.extent = sub.ext){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.crop <- raster::crop(x[[1]],y = raster::extent(crop.extent))  
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  # convert to a data frame and add in a half month at average of Jan+Dec at either end
  colnames(x.cv)[3:ncol(x.cv)] <- 1:12
  x.df <- x.cv[complete.cases(x.cv),]
  x.df$`0.5` <- x.df$`12.5` <- (x.df$`1` + x.df$`12`)/2
  # gather everything together and compute averages per month across the latitude breaks
  ER.df <- x.df %>%
    tidyr::gather("month","ER",3:ncol(x.df), convert=TRUE) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(ER.mean = mean(ER))

  # rolling average for fill colour
  ER.df2 <- ER.df %>% group_by(strata) %>% arrange(month) %>% mutate(Roll.mean = roll_mean(ER.mean, n=2, fill=0))

  # plot
  ER.ridge <- ggplot(data= ER.df2,
         aes(x= month,y= strata,height = ER.mean, group = strata, fill = Roll.mean)) +
    geom_density_ridges_gradient(color="#0000000F", stat = "identity", scale = scale) +
    scale_fill_gradient(low= "yellow", high = "red4",
                        limits = c(0,max(ER.df$ER.mean)),
                        name = "Mean \nEbola \nRisk") +
    scale_x_continuous(breaks = 1:12, labels=month.abb[1:12],
                       expand = c(0,0))+
    scale_y_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour="grey96")
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

