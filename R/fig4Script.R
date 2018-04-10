##########################
#### Figure 4         ####
##########################
source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster)
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
mylog = scales::trans_new('mylog',
                          transform=function(x) { log(x+1e-5) },
                          inverse=function(x) { exp(x)-1e-5},
                          breaks = scales::log_breaks(base=exp(1)),
                          domain=c(1e-20, 100))



afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")
rf.poly <- rasterToPolygons(raster(file.path(data.source, "cropMask.tif")),
                            fun = function(x){x==1}, dissolve = T)
bkg <- theme(
  plot.title = element_text(hjust = 0.5),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())

scaleFUN <- function(x)sprintf("%.4f", round(x, digits = 4)) 
                                
ERgplot <- function(x, source.path = data.source, afr = afr.poly, rf = rf.poly){
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Risk")
  
  g.plot <- ggplot(sum.df) +
    aes(x=long, y=lat) +
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 alpha = .2) +
    
    #fill data values
    geom_raster(aes(fill = Risk), interpolate = T)+
    scale_fill_gradientn(name = "Average \nRisk", trans = scales::log_trans(), na.value=terrain.colors(10)[1],
                         colors = terrain.colors(10), limits = c(1e-4, 12)) +
   
    #add area modeled
    geom_polygon(data = fortify(rf),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15)) +
    scale_y_continuous(expand = c(0,0), lables = scaleFUN) +
    theme_bw()+
    theme(axis.title = element_blank())
    
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15))
  out <- g.plot + theme_bw() +bkg
}


risk.plot <- ERgplot(hum.mean)
risk.plot

ggsave("figures/fig4_A.pdf",
       risk.plot,
       device = "pdf",
       width = 5,
       height = 5,
       units = "in")

#### Pannel 2 Facetted monthly ####
#nullFacet theme
nullFacet <- theme_minimal() + theme(
  #remove all axis info... it's for the birds
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  #get them in close (real close)
  panel.spacing.x = unit(-1, "lines"),
  panel.spacing.y = unit(-.5, "lines"),
  #hack off names and put them in the Atlantic somewhere
  strip.text = element_blank()
)

facetRisk <- function(x, source.path = data.source, afr= afr.poly, rf = rf.poly){
  ## results dataframe
  res <- data.frame(rasterToPoints(x[[1]]))
  colnames(res) <- c("long","lat","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  res.long <- tidyr::gather(data = res, key = "Month", value = "Risk", Jan:Dec, factor_key = T)

  monthly.plot <- ggplot(res.long) +
    ### Monthly rasters
    aes(x=long, y=lat) +
    geom_raster(aes(fill = Risk), interpolate = T)+
    scale_fill_gradientn(name = "Risk", trans = scales::log_trans(), na.value=terrain.colors(10)[1],
                         colors = terrain.colors(10), limits = c(1e-4, 12)) +
    #add area modeled
    geom_polygon(data = fortify(rf),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 fill = NA,
                 alpha = .1) +
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15)) +
    facet_wrap(~ Month, ncol = 3)+
    nullFacet 
  
  return(monthly.plot)
}
z <- facetRisk(hum.mean)

ggsave("figures/fig4_B.png",
       z,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
