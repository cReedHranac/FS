##################
#### SI Plots ####
##################

source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

#### functions ####
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

birthForce1 <- function(x, afr = afr.poly, save = F, crop.extent = Afr.ext, device = NULL, ....){
  ###Function for creating facted birthforce maps
  ##gen dataframe
  stk.crop <- crop(x[[1]], crop.extent)
  stk.df <- data.frame(rasterToPoints(stk.crop))
  #create better names
  base <- substring(names(x[[1]])[1], 4, 6) 
  i <- 1:12
  j <- c(i[2:12], i[1])
  colnames(stk.df) <- c("long", "lat", paste0(base,i,"_",j))
  res.long <- tidyr::gather(data = stk.df, key = "window", value = "BF", -c(long, lat), factor_key = T)
  
  ##plot
  bf1.plot <- ggplot(res.long, aes(x = long, y = lat)) +
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = NA,
                 fill = 'black',
                 alpha = .2) +
    
    #fill Raster values
    geom_raster(aes(fill = BF), interpolate = T)+
    
    #Colors
    scale_fill_gradientn(colors = c("#5e3c99", "#b2abd2","#ffffff", "#fdb863","#e66101"), 
                         limits = c(0,max(res.long$BF)),
                         name = "Number \nBirthing")+
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "grey20",
                 alpha = .20) +
    
    # limit coordinates
    coord_fixed(xlim = crop.extent[1:2],ylim = crop.extent[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    #Extras
    theme_bw() + 
    theme( axis.title = element_blank()) +
    scale_y_continuous(expand = c(0,0)) +
    
    #Facet  
    facet_wrap(~ window, ncol= 3)
  
  if(save == T){
    ggsave(filename = file.path("figures/", paste0(base,"BF_SI", ".", device)),
           bf1.plot, device = device)
  }
  
  return(bf1.plot)
}

#### Monthy Force of Breeding ####
ptr.sum <- sumGen(model.string = "ptr.dbl.imp")
mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")

Afr.ext <- c(-18, 49, -36, 16)
crop.extent <- Afr.ext

ptr.bf1 <- birthForce1(ptr.sum, save = T, device = "png")
mic.bf1 <- birthForce1(mic.sum, save = T, device = "png")
mol.bf1 <- birthForce1(mol.sum, save = T, device = "png")