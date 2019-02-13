##################
#### SI Plots ####
##################

source("R/helperFunctions.R")
library(tidyverse); library(gtools)
library(RcppRoll);library(ggridges)

# default theme for ggplot
theme_set(theme_bw(base_size = 9, base_family = "serif"))

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

better.names <- function(x, model.name){
  ### function for impoving names accociated with items retrieved from SpatHandler
  base <- sapply(strsplit(model.name, "_"), tail, 1)
  i <- 1:12
  j <- c(i[12],i[1:11])
  names(x) <- paste0(base, "_", j, "_", i)
  return(x)
}

spatHandler.Proj <- function(model.name, mod.dir){
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
  a <- stack(f.list)
  stk <- better.names(a, model.name )
  m.stk <- mean(stk)
  out.l <- list(stk,m.stk)
  return(out.l)
}

ERgplot <- function(x, source.path = data.source, afr = afr.poly, rf = rf.poly){
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Risk")
  
  #create african continent background
  g.plot <- ggplot(sum.df, aes(x=long, y=lat)) +
    
    # background
    geom_polygon(data = fortify(afr),
                 aes(group = group), 
                 colour = NA,
                 fill = 'black',
                 alpha = .2) +
    
    # fill data values
    geom_raster(aes(fill = Risk), interpolate=TRUE) +
    
    # color scale
    scale_fill_gradient2(trans = scales::log_trans(),
                         low = "yellow", mid = "red4",
                         #limits = c(1e-20, 2),
                         limits = c(1e-4, 12),
                         na.value = "yellow", 
                         name = "Mean \nEbola \nRisk") +
    
    # country borders
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "black",
                 fill = NA,
                 size=0.1) +
    
    
    # limit coordinates
    coord_fixed(xlim = Africa.ext[1:2], ylim = Africa.ext[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    # theme
 #   theme_bw() +
    theme(axis.title = element_blank())
  
  return(g.plot)
}

ERridge <- function(x, n.bin, scale = 5, crop.extent = Africa.ext){
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
  ER.df2 <- ER.df %>% 
    group_by(strata) %>%
    arrange(month) %>%
    mutate(Roll.mean = roll_mean(ER.mean, n=2, fill=0))
  
  # plot
  ER.ridge <- ggplot(data= ER.df2,
                     aes(x= month,y= strata,height = ER.mean, group = strata, fill = Roll.mean)) +
    geom_density_ridges_gradient(color="#0000000F", stat = "identity", scale = scale) +
    scale_fill_gradient(low= "yellow", high = "red4",
                        limits = c(0,max(ER.df$ER.mean)),
                        name = "Mean \nEbola \nRisk") +
    scale_x_continuous(breaks = 1:12, labels=substring(month.abb, 1, 1),
                       expand = c(0,0))+
    scale_y_discrete(expand=c(0,0)) +
  #  theme_bw() +
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

birthForce1 <- function(x, c.string,  afr = afr.poly, save = F, crop.extent = Africa.ext, device.out = NULL, dir.out, ....){
  ###Function for creating facted birthforce maps
  ##gen dataframe
  stk.crop <- crop(x[[1]], crop.extent)
  stk.df <- data.frame(rasterToPoints(stk.crop))
  #create better names
  base <- substring(names(x[[1]])[1], 4, 6) 
  i <- 1:12
  j <- c(i[2:12], i[1])
  colnames(stk.df) <- c("long", "lat", paste0(base,i,"_",j))
  res.long <- tidyr::gather(data = stk.df, key = "window", value = "BF", -c(long, lat)) %>%
    tidyr::extract(window, into=c("first", "second"), regex='([0-9]+)_([0-9]+)', convert=TRUE) %>%
    mutate(window = paste(month.abb[first], month.abb[second], sep='-')) %>%
    mutate(window = factor(window, levels = paste(month.abb, month.abb[c(2:12,1)], sep='-')))
  
  ##plot
  bf1.plot <- ggplot(res.long, aes(x = long, y = lat)) +
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "#212121",
                 fill = NA) +
    
    #fill Raster values
    geom_raster(aes(fill = BF), interpolate = T)+
    
    #Colors
    scale_fill_gradientn(colors = c.string, 
                         limits = c(0,max(res.long$BF)),
                         name = "Number \nBirthing")+
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "#212121",
                 fill = NA) +
    aes(x = long, y = lat)+
    
    # limit coordinates
    coord_fixed(xlim = crop.extent[1:2],ylim = crop.extent[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    #Extras
#    theme_bw() + 
    theme( axis.title = element_blank(),
           strip.text.x = element_text(margin = margin(1, 0, 1, 0, "pt"))) +
    
    #Facet  
    facet_wrap(~ window, ncol= 3)
  
  if(save == T){
    if(device.out == "pdf"){
      dev.ext <- cairo_pdf
    } else if (device.out =="eps"){
      dev.ext <- cairo_ps
      } else {
      dev.ext <- device.out
    }
    ggsave(filename = file.path(dir.out, paste0(base,"BF_SI", ".", device.out)),
           bf1.plot, width = 6, height = 6.5, units = "in", device = dev.ext)
  }
  
  return(bf1.plot)
}

riskForce1 <- function(x, afr = afr.poly, save = F, crop.extent = Africa.ext, device.out = NULL, dir.out, file.out, ....){
  ###Function for monthly facet plotting of Human and Animal Risk
  stk.crop <- crop(x[[1]], crop.extent)
  stk.df <- data.frame(rasterToPoints(stk.crop))
  ##creating better names
  base <- substring(names(x[[1]]), 1, 3)
  i <- 1:12
  j <- c(i[12],i[1:11])
  colnames(stk.df) <- c("long", "lat", paste0(base,"_", j,"_", i))
  ##modify structure for our purposes
  res.long <- tidyr::gather(data = stk.df, key = "window", value = "Risk", -c(long, lat)) %>%
    tidyr::extract(window, into=c("first", "second"), regex='([0-9]+)_([0-9]+)', convert=TRUE) %>%
    mutate(window = paste(month.abb[first], month.abb[second], sep='-')) %>%
    mutate(window = factor(window, levels = paste(month.abb, month.abb[c(2:12,1)], sep='-')))

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
                         labels = scales::number_format(accuracy=0.001),
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
 #   theme_bw() + 
    theme( axis.title = element_blank(),
           strip.text.x = element_text(margin = margin(1, 0, 1, 0, "pt"))) +

    #Facet  
    facet_wrap(~ window, ncol= 3)
  
  if(save == T){
    if(device.out == "pdf"){
      dev.ext <- cairo_pdf
    } else if (device.out =="eps"){
      dev.ext <- cairo_ps
    } else {
      dev.ext <- device.out
    }
    if (missing(file.out))
      file.out = paste0(base,"Risk_SI", ".", device.out)
    ggsave(filename = file.path(dir.out, file.out),
           risk1.plot, width = 6, height = 6.5, units = "in", device = dev.ext)
  }
  
  return(risk1.plot)
  
}

#### Data for plotting ####
## DRC health districs
library(rgdal)
afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")

Africa.ext <- c(-18, 47, -36, 16)

drc.hd <- readOGR(file.path(data.source, "zone_ste_puc"),
                  "Zone_Ste_Puc")
## proj4 I use
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
drc <- spTransform(drc.hd, CRS(wgs))
rm(drc.hd)


#### Animal Top outbreak ####
ann.mean <- spatHandler.Proj(model.name = "a_Prb",mod.dir =  dOut.an)
ann.risk.plot <- ERgplot(ann.mean)
n.ridge <- 96  # Number of ridgelines. Mess with scale below as well.
w.ridge <- 0.5 # Width of ridge plot compared to map
ann.ridge <- ERridge(ann.mean, n.bin = n.ridge, scale = 2, crop.extent = Africa.ext )
# and fixup the aspect ratio of ridge to that of map
width_height <- diff(Africa.ext)[c(1,3)]
aspect_map <- width_height[1] / width_height[2]

# now aspect ratio of ridge plot (this assumes top bin isn't too large)
aspect_ridge <- 12 / n.ridge
aspect_match <- aspect_ridge / aspect_map / w.ridge

ANrisk.grob <- ggplotGrob(ann.risk.plot + guides(fill='none') +
                            theme(plot.margin=unit(rep(0,4), "cm")))
ANridge.grob <- ggplotGrob(ann.ridge + guides(fill='none') +
                             theme(plot.margin=unit(rep(0,4), "cm")) +
                             coord_fixed(ratio = aspect_match))

# find the height and left axis width of the risk grob
index <- ANrisk.grob$layout$t[ANrisk.grob$layout$name == 'panel']
plot_height <- ANrisk.grob$heights[index]
index <- ANrisk.grob$layout$l[ANrisk.grob$layout$name == 'axis-l']
axis_width <- ANrisk.grob$widths[index]
index <- ANrisk.grob$layout$l[ANrisk.grob$layout$name == 'axis-r']
ANrisk.grob$widths[index] <- unit(0.5, 'cm')

# set the height the same, the left axis width, and the right
# axis to the same as the left axis for the map, scaled so that
# the widths are maintained (when laying out with grid.arrange)
# the widths are applied to the full objects, so we want everything
# to be a scaled version of the other (i.e. axis placement/widths)
# etc scaled perfectly
index <- ANridge.grob$layout$t[ANridge.grob$layout$name == 'panel']
ANridge.grob$heights[index] <- unit(as.numeric(plot_height)/w.ridge, 'null')
index <- ANridge.grob$layout$l[ANridge.grob$layout$name == 'axis-l']
right_axis_width <- unit(
  grid::convertWidth(axis_width, unitTo = 'cm', valueOnly = TRUE) * w.ridge,
  "cm")
ANridge.grob$widths[index] <- right_axis_width
index <- ANridge.grob$layout$l[ANridge.grob$layout$name == 'axis-r']
ANridge.grob$widths[index] <- unit(0.5 * w.ridge, 'cm')

library(gridExtra)
# arrange them side by side
an.risk <- grid.arrange(ANrisk.grob,
                        ANridge.grob,
                        ncol=2, widths=c(1,w.ridge))

ggsave(filename = file.path(fig.si,"AnRisk_TopModel.pdf"),
       an.risk,
       device = cairo_pdf,
       width=8,
       height=4.3,
       units = "in",
       dpi = 300)




library(raster)
#### Monthy Force of Breeding ####
ptr.sum <- sumGen(model.string = "ptr.dbl.imp")
## fix ptr namess
c.names <- names(ptr.sum[[1]])
fx.names <- gsub("ptr", "afb", c.names)
names(ptr.sum[[1]]) <- fx.names

mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")

c1 <- colorRampPalette(c("#CFD8DC", "green4"))
c2 <- colorRampPalette(c("#CFD8DC", "dodgerblue2"))
c3 <- colorRampPalette(c("#CFD8DC", "darkorange2"))
library(ggridges)
ptr.bf1 <- birthForce1(x = ptr.sum,
                       c.string = c1(5) ,
                       save = T,
                       device.out = "pdf",
                       dir.out = fig.si)
mol.bf1 <- birthForce1(mol.sum,
                       c.string = c2(5),
                       save = T,
                       device.out = "pdf",
                       dir.out = fig.si)
mic.bf1 <- birthForce1(mic.sum,
                       c.string = c3(5),
                       save = T,
                       device.out = "pdf",
                       dir.out = fig.si)

#### Monthly Human and Animal Risk Maps
hum.sum <- spatHandler.Proj(model.name = "h_Prb",mod.dir =  dOut.1)

riskForce1(hum.sum,
           save = T,
           device.out = "pdf",
           file.out = "HumRisk_SI.pdf",
           dir.out = fig.si)

ann.sum <- spatHandler.Proj("a_Prb", mod.dir = dOut.an)
riskForce1(ann.sum,
           save = T,
           device.out = "pdf",
           file.out = "AnnRisk_SI.pdf",
           dir.out = fig.si)
