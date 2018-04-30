###############################
#### Figure 3              ####
###############################
source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra)
library(raster); library(rgdal); library(tidyverse)

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
ptr.sum <- sumGen(model.string = "ptr.dbl.imp")
mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")


afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")
rf.poly <- rasterToPolygons(raster(file.path(data.source, "cropMask.tif")),
                            fun = function(x){x==1}, dissolve = T)
bkg <- theme(
  plot.title = element_text(hjust = 0.5),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())

BFgplot <- function(x, afr = afr.poly, rf = rf.poly, themed = bkg){
  #### Set Up ####
  afr.poly <- afr
  rf.poly <- rf
  
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Number")
  
  g.plot <- ggplot(sum.df) +
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 alpha = .20) +
    aes(x=long, y=lat) +
    
    #Colors
    scale_fill_gradient(low = "#f4eade", high = "#2988bc",
                        limits = c(0,max(sum.df$Number)),
                        name = "Number \nBirthing")+
    #Raster
    geom_raster(aes(fill = Number), interpolate = T)+
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 fill = NA,
                 alpha = .2) +
    
    #add area modeled
    geom_polygon(data = fortify(rf.poly),
                 aes(long, lat, group = group),
                 colour = "white", 
                 fill = NA) +
    
    #Extras
    coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15)) + 
    theme_bw() + 
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = c(.2,.4),
           legend.background = element_blank())+
  scale_y_continuous(expand = c(0,0))
  
}

ptr.BF <- BFgplot(x = ptr.sum)
mol.BF <- BFgplot(x = mol.sum)
mic.BF <- BFgplot(x = mic.sum)

ggsave("figures/fig3_A.png",
       ptr.BF,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
ggsave("figures/fig3_B.png",
       mol.BF,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
ggsave("figures/fig3_C.png",
       mic.BF,
       device = "png",
       width = 5,
       height = 5,
       units = "in")


#### Pannel 2 (Right) ####
# install_github("cran/ggridges")
library(ggridges); library(readr); library(tidyr)
sub.ext <- c(-18, 49, -36, 15) #extent subset like that of the other map figures
BFridge <- function(x, n.bin, crop.extent = sub.ext){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.crop <- crop(x[[1]], crop.extent)
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  colnames(x.cv)[3:ncol(x.cv)] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  x.df <- x.cv[complete.cases(x.cv),]
  x.df$Jan2 <- x.df$Jan
  bf.df <- x.df %>%
    tidyr::gather("month","BF",3:ncol(x.df)) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(bf.mean = mean(BF)) 
  
  bf.df$month <-  factor(bf.df$month,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec", "Jan2"))
  
  
  bf.ridge <- ggplot(data= bf.df, 
                     aes(x= month,y= strata,height = bf.mean, group = strata, fill = bf.mean))+
    geom_density_ridges_gradient(stat = "identity", scale = 3, alpha = .5, aes()) +
    scale_fill_gradient(low = "#f4eade", high = "#2988bc",
                        limits = c(0,max(bf.df$bf.mean)),
                        name = "Mean \nBirth \nForce") +
    scale_x_discrete(label = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec", "Jan"),
                     expand = c(0,0))+
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
       )
  
  return(bf.ridge )
}

r.ptr <- BFridge(x = ptr.sum, n.bin = 40, crop.extent = sub.ext)
r.mic <- BFridge(mic.sum, 40)
r.mol <- BFridge(mol.sum, 40)

ggsave("figures/fig3_D.png",
       r.ptr,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
ggsave("figures/fig3_E.png",
       r.mol,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
ggsave("figures/fig3_F.png",
       r.mic,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
fig3.complete <- grid.arrange(ptr.BF, r.ptr, 
             mic.BF, r.mic,
             mol.BF, r.mol, 
             layout_matrix= rbind(c(1,2),
                                  c(3,4),
                                  c(5,6)))
ggsave("figures/Fig3Complete.png",
      fig3.complete,
      device = "png", 
      width = 7.5,
      height = 7.5,
      units = "in")
