##################
#### SI Plots ####
##################

source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

#### functions ####
sumGen <- function(model.string){
  ## function for loading rasters and producing an averaged product based on the USED FOR BirthFORCE ITEMS
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

better.names <- function(x){
  ### function for impoving names accociated with items retrieved from SpatHandler
  base <- substring(names(x[[1]]), 1, 3)
  i <- 1:12
  j <- c(i[12],i[1:11])
  names(x) <- paste0(base, "_", j, "_", i)
  return(x)
}

birthForce1 <- function(x, afr = afr.poly, save = F, crop.extent = Afr.ext, device.out = NULL, ....){
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
                 fill = "white",
                 alpha = .2) +
    
    #fill Raster values
    geom_raster(aes(fill = BF), interpolate = T)+
    
    #Colors
    scale_fill_gradientn(colors = c("#CFD8DC","#BDBDBD", "#00BCD4", "#0097A7"), 
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
    if(device.out == "pdf"){
      dev.ext <- cairo_pdf
    }else{
      dev.ext <- device.out
    }
    ggsave(filename = file.path("figures/", paste0(base,"BF_SI", ".", device.out)),
           bf1.plot, device = dev.ext)
  }
  
  return(bf1.plot)
}

riskForce1 <- function(x, afr = afr.poly, save = F, crop.extent = Afr.ext, device.out = NULL, ....){
  ###Function for monthly facet plotting of Human and Animal Risk
  stk.crop <- crop(x[[1]], crop.extent)
  stk.df <- data.frame(rasterToPoints(stk.crop))
  ##creating better names
  base <- substring(names(x[[1]]), 1, 3)
  i <- 1:12
  j <- c(i[12],i[1:11])
  colnames(stk.df) <- c("long", "lat", paste0(base,"_", j,"_", i))
  ##modify structure for our purposes
  res.long <- tidyr::gather(data = stk.df, key = "window", value = "Risk", -c(long, lat), factor_key = T)
  
  ##plot
  risk1.plot <- ggplot(res.long, aes(x = long, y = lat)) +
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = NA,
                 fill = 'black',
                 alpha = .2) +
    
    #fill Raster values
    geom_raster(aes(fill = Risk), interpolate = T)+
    
    #Colors
    scale_fill_gradient2(trans = scales::log_trans(),
                         low = "yellow", mid = "red4",
                         limits = c(1e-4, 12),
                         na.value = "yellow", 
                         name = "Ebola \nSpillover \nRisk") +
    
    
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
    if(device.out == "pdf"){
      dev.ext <- cairo_pdf
    }else{
      dev.ext <- device.out
    }
    ggsave(filename = file.path("figures/", paste0(base,"Risk_SI", ".", device.out)),
           risk1.plot, device = dev.ext)
  }
  
  return(risk1.plot)
  
}

#### Data for plotting ####
afr.poly <- readOGR(dsn = path.expand(file.path(data.source, "Africa")),
                    layer = "AfricanCountires")
Afr.ext <- c(-18, 49, -36, 16)


#### Monthy Force of Breeding ####
ptr.sum <- sumGen(model.string = "ptr.dbl.imp")
mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")



ptr.bf1 <- birthForce1(ptr.sum, save = T, device.out = "pdf")
mic.bf1 <- birthForce1(mic.sum, save = T, device.out = "pdf")
mol.bf1 <- birthForce1(mol.sum, save = T, device.out = "pdf")

#### Monthly Human and Animal Risk Maps
hum.sum <- spatHandler("humNoAn", "SpGLMRes_F")
riskForce1(hum.sum, save = T, device.out = "pdf")

ann.sum <- spatHandler("ann", "SpGLMRes_F")
riskForce1(ann.sum, save = T, device.out = "pdf")
