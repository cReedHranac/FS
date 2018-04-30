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
    scale_fill_gradient(name = "Average \nRisk", trans = scales::log_trans(), na.value = "yellow",
                         low = "yellow", high = "red4", limits = c(1e-4, 12)) +
   
    #add area modeled
    geom_polygon(data = fortify(rf),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15)) +
    scale_y_continuous(expand = c(0,0) )#, lables = scaleFUN) +
    theme_bw()+
    theme(axis.title = element_blank())
    
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15))
  out <- g.plot + theme_bw() +bkg
}


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


#### Idea for zoom pannels ####
insert.ext <- c(-15, 40.5,-5,10)

highRiskZone <- risk.plot + coord_fixed(xlim = insert.ext[1:2],ylim = insert.ext[3:4])

library(ggridges)

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
k <- ERridge(x = hum.mean, n.bin = 100, crop.extent = insert.ext)
k

ggsave("figures/fig4_B.png",
       highRiskZone,
       device = "png",
       width = 5,
       height = 5,
       units = "in")


ggsave("figures/fig4_c.png",
       k,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
grid.arrange(risk.plot, highRiskZone, k, 
             layout_matrix = rbind(c(1,1,2,2,3,3),
                                   c(1,1,2,2,3,3)))




#### Alternative human models for ridges ####
hum.bat <- spatHandler("humBat")
batRidge <- ERridge(hum.bat, n.bin = 100, insert.ext)

hum.NoAn <- spatHandler("humNoAn") # Go with this one for ridges
NoAn.risk <- ERgplot(hum.NoAn)
NoAnRidge <- ERridge(hum.NoAn, n.bin = 100, insert.ext)

#### More Items ####
hum.mean <- spatHandler("hum") #This for averages 
risk.plot <- ERgplot(hum.mean)


ggsave("figures/fig4_A.png",
       risk.plot,
       device = "png",
       width = 5,
       height = 5,
       units = "in")

Afr.ext <- c(-18, 49, -36, 15)
afr.ridge <- ERridge(hum.NoAn, n.bin = 50, crop.extent = Afr.ext )

ggsave("figures/fig4_B.png",
       afr.ridge,
       device = "png",
       width = 5,
       height = 5,
       units = "in")


centeral.africa <- c(8,35,-5,6)
CentralZone <- risk.plot + coord_fixed(xlim = centeral.africa[1:2],ylim = centeral.africa[3:4])
central.ridge <- ERridge(hum.NoAn, n.bin = 50, centeral.africa)
ggsave("figures/fig4_C.png",
       CentralZone,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
ggsave("figures/fig4_D.png",
       central.ridge,
       device = "png",
       width = 5,
       height = 5,
       units = "in")


west.africa <- c(-14,5, 4, 12)
westernZone <- risk.plot + coord_fixed(xlim = west.africa[1:2],ylim = west.africa[3:4])
western.ridge <- ERridge(hum.NoAn, n.bin = 100, west.africa)
ggsave("figures/fig4_E.png",
       westernZone,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
ggsave("figures/fig4_F.png",
       western.ridge,
       device = "png",
       width = 5,
       height = 5,
       units = "in")


fig4.complete <- grid.arrange(risk.plot, afr.ridge,
             CentralZone, central.ridge,
             westernZone, western.ridge,
            layout_matrix = rbind(c(1,3,5),
                                  c(1,3,5),
                                  c(2,4,6)))
ggsave("figures/fig4_Complete.png",
       fig4.complete,
       device = "png",
       width = 7,
       height = 7,
       units = "in")
