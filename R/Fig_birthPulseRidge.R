###############################
#### Figure 3              ####
###############################
source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

theme_set(theme_bw(base_family="serif"))

# Africa extent to use. This is reasonably tight around the data
Africa.ext <- c(-18, 47, -36, 16)

#### Functions ####

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
BFgplot <- function(x, afr = afr.poly, c.string){
  ## dataframe for plotting
  sum.df <- data.frame(rasterToPoints(x[[2]]))
  colnames(sum.df) <- c("long","lat","Number")
  
  g.plot <- ggplot(sum.df, aes(x=long, y=lat)) +
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = NA,
                 fill = "#FFFFFF",
                 alpha = .2) +
    
    #fill Raster values
    geom_raster(aes(fill = Number), interpolate = T)+
    
    #Colors
    scale_fill_gradientn(colors = c.string, 
                        limits = c(0,max(sum.df$Number)),
                        name = "Mean \nNumber \nBirthing")+

    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "#212121",
                 fill = NA) +
    aes(x=long, y=lat) +
    
    # limit coordinates
    coord_fixed(xlim = Africa.ext[1:2], ylim = Africa.ext[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    #Extras
    theme( axis.title = element_blank(),
           legend.position = c(0.2,0.3)) +
    scale_y_continuous(expand = c(0,0))
  
  return(g.plot)
}
BFridge <- function(x, n.bin, crop.extent = Africa.ext, scale, c.string){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.crop <- crop(x[[1]], crop.extent)
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  # convert to a data frame and add in a half month at average of Jan+Dec at either end
  colnames(x.cv)[3:ncol(x.cv)] <- 1:12
  x.df <- x.cv[complete.cases(x.cv),]
  x.df$`0.5` <- x.df$`12.5` <- (x.df$`1` + x.df$`12`)/2
  # gather everything together and compute averages per month across the latitude breaks
  bf.df <- x.df %>%
    tidyr::gather("month","BF",3:ncol(x.df), convert=TRUE) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(bf.mean = mean(BF)) 
  
  # rolling average for fill colour
  bf.df2 <- bf.df %>% 
    group_by(strata) %>%
    arrange(month) %>%
    mutate(Roll.mean = roll_mean(bf.mean, n=2, fill=0))
  
  
  #Plot
  bf.ridge <- ggplot(data= bf.df2, 
                     aes(x= month,y= strata,height = bf.mean, group = strata, fill = Roll.mean)) +
    geom_density_ridges_gradient(color="#0000000F", stat = "identity", scale=scale) +
    scale_fill_gradientn(colors = c.string,
                        limits = c(0,max(bf.df2$bf.mean)),
                        name = "Mean \nBirth \nForce") +
    scale_x_continuous(breaks = 1:12, labels=substring(month.abb, 1, 1),
                       expand = c(0,0))+
    scale_y_discrete(expand=c(0,0)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour="grey96")
    )
  
  return(bf.ridge )
}

#### Left Pannels ####
afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")

ptr.sum <- sumGen(model.string = "ptr.dbl.imp")
mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")

## Create color ramp palets to match with figure 2
c1 <- colorRampPalette(c("#CFD8DC", "green4"))
c2 <- colorRampPalette(c("#CFD8DC", "dodgerblue2"))
c3 <- colorRampPalette(c("#CFD8DC", "darkorange2"))

ptr.BF <- BFgplot(x = ptr.sum, c.string = c1(5))
mol.BF <- BFgplot(x = mol.sum, c.string = c2(5))
mic.BF <- BFgplot(x = mic.sum, c.string = c3(5))

#### Pannel 2 (Right) ####
n.ridge <- 48  # Number of ridgelines. Mess with scale below as well.
w.ridge <- .5 # Width of ridge plot compared to map

r.ptr <- BFridge(x = ptr.sum,
                 n.bin = n.ridge,
                 scale = 2,
                 crop.extent = Africa.ext,
                 c.string = c1(5))
r.mol <- BFridge(mol.sum,
                 n.bin = n.ridge,scale = 2,
                 crop.extent = Africa.ext, 
                 c.string = c2(5))
r.mic <- BFridge(mic.sum,
                 n.bin = n.ridge,scale = 2,
                 crop.extent = Africa.ext, 
                 c.string = c3(5))

# ggsave("figures/fig3_D.png",
#        r.ptr,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in")
# ggsave("figures/fig3_E.png",
#        r.mol,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in")
# ggsave("figures/fig3_F.png",
#        r.mic,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in")

# Arrange the map and ridge on the same plot
# first work out aspect ratio of map
width_height <- diff(Africa.ext)[c(1,3)]
aspect_map <- width_height[1] / width_height[2]

# now aspect ratio of ridge plot (this assumes top bin isn't too large)
aspect_ridge <- 12 / n.ridge

# and fixup the aspect ratio of ridge to that of map
aspect_match <- aspect_ridge / aspect_map / w.ridge

map.list <- list(ptr.BF, mol.BF, mic.BF)
ridge.list <- list(r.ptr, r.mol, r.mic)

map.grob <- lapply(map.list, function(x){ggplotGrob(x + guides(fill='none') +
                                                     theme(plot.margin=unit(rep(0,4), "cm")))})
ridge.grob <- lapply(ridge.list, function(x){ggplotGrob(x + guides(fill='none') +
                                                          theme(plot.margin=unit(rep(0,4), "cm"))+
                                                          coord_fixed(ratio = aspect_match))})

# find the height and left axis width of the risk grob

risk.margin.reset <- function(x){
  index <- x$layout$t[x$layout$name == 'panel']
  plot_height <- x$heights[index]
  index <- x$layout$l[x$layout$name == 'axis-l']
  axis_width <- x$widths[index]
  index <- x$layout$l[x$layout$name == 'axis-r']
  x$widths[index] <- unit(0.5, 'cm')
  return(x)
}

risk.list <- lapply(map.grob, risk.margin.reset)

# set the height the same, the left axis width, and the right
# axis to the same as the left axis for the map, scaled so that
# the widths are maintained (when laying out with grid.arrange)
# the widths are applied to the full objects, so we want everything
# to be a scaled version of the other (i.e. axis placement/widths)
# etc scaled perfectly

# Grab the correct axis width/plot heights from the first map panel
x <- map.grob[[1]]
index <- x$layout$t[x$layout$name == 'panel']
plot_height <- x$heights[index]
index <- x$layout$l[x$layout$name == 'axis-l']
axis_width <- x$widths[index]
index <- x$layout$l[x$layout$name == 'axis-r']
x$widths[index] <- unit(0.5, 'cm')

ridge.margit.reset <- function(x, plot_height, axis_width) {
  index <- x$layout$t[x$layout$name == 'panel']
  x$heights[index] <- unit(as.numeric(plot_height)/w.ridge, 'null')
  index <- x$layout$l[x$layout$name == 'axis-l']
  right_axis_width <- unit(
    grid::convertWidth(axis_width, unitTo = 'cm', valueOnly = TRUE) * w.ridge,
    "cm")
  x$widths[index] <- right_axis_width
  index <- x$layout$l[x$layout$name == 'axis-r']
  x$widths[index] <- unit(0.5 * w.ridge, 'cm')
  # x$respect <- TRUE
  return(x)
}

ridge.list <- lapply(ridge.grob, ridge.margit.reset, plot_height=plot_height, axis_width=axis_width)

master.list <- list(risk.list[[1]], ridge.list[[1]],
                    risk.list[[2]], ridge.list[[2]],
                    risk.list[[3]], ridge.list[[3]])

#making individual ones for presentations
# afrP <- grid.arrange(grobs = master.list[1:2], 
#                      ncol = 2, widths=c(1,w.ridge)) 
# micP <- grid.arrange(grobs = master.list[5:6], 
#                      ncol = 2, widths=c(1,w.ridge)) 
# molp <- grid.arrange(grobs = master.list[3:4], 
#                              ncol = 2, widths=c(1,w.ridge)) 
# # ggsave("figures/afbPulse.pdf",
#        afrP,
#        device = cairo_pdf,
#        width=6,
#        height=3.2,
#        units = "in",
#        dpi = 300)
# ggsave("figures/molPulse.pdf",
#        molp,
#        device = cairo_pdf,
#        width=6,
#        height=3.2,
#        units = "in",
#        dpi = 300)
# ggsave("figures/micPulse.pdf",
#        micP,
#        device = cairo_pdf,
#        width=6,
#        height=3.2,
#        units = "in",
#        dpi = 300)




f3 <- grid.arrange(grobs = master.list,
             ncol = 2, widths=c(1,w.ridge))
ggsave(file.path(fig.pub,"birthForceComplete.pdf"),
       f3,
       device = cairo_pdf,
       width=6,
       height=9.6,
       units = "in",
       dpi = 300)



#### Create the same thing for pure probability ####
BFgplotP <- function(x, afr = afr.poly, c.string){
  ## dataframe for plotting
  x.prob <- calc(x[[1]], function(x){ifelse(is.na(x),NA,x/1000)})
  sum.df <- data.frame(rasterToPoints(x.prob))
  colnames(sum.df) <- c("long","lat","Number")
  
  g.plot <- ggplot(sum.df, aes(x=long, y=lat)) +
    
    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = NA,
                 fill = "#FFFFFF",
                 alpha = .2) +
    
    #fill Raster values
    geom_raster(aes(fill = Number), interpolate = T)+
    
    #Colors
    scale_fill_gradientn(colors = c.string, 
                         limits = c(0,1),
                         name = "Average\nBirthing\nProbability")+
    
    #create african continent background
    geom_polygon(data = fortify(afr),
                 aes(long, lat, group = group), 
                 colour = "#212121",
                 fill = NA) +
    aes(x=long, y=lat) +
    
    # limit coordinates
    coord_fixed(xlim = Africa.ext[1:2], ylim = Africa.ext[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    #Extras
    theme( axis.title = element_blank(),
           legend.position = c(0.2,0.3)) +
    scale_y_continuous(expand = c(0,0))
  
  return(g.plot)
}
BFridgeP <- function(x, n.bin, crop.extent = Africa.ext, scale, c.string){
  ## Function for creating ridgeline density plots of the breeding force
  ## used on objects creaded from sumGen (since it loads rasterlayers as well)
  x.fix <- calc(x[[2]], function(x){ifelse(is.na(x),NA,x/1000)})
  x.crop <- crop(x.fix, crop.extent)
  x.cv <- as.data.frame(rasterToPoints(x.crop))
  # convert to a data frame and add in a half month at average of Jan+Dec at either end
  colnames(x.cv)[3:ncol(x.cv)] <- 1:12
  x.df <- x.cv[complete.cases(x.cv),]
  x.df$`0.5` <- x.df$`12.5` <- (x.df$`1` + x.df$`12`)/2
  # gather everything together and compute averages per month across the latitude breaks
  bf.df <- x.df %>%
    tidyr::gather("month","BF",3:ncol(x.df), convert=TRUE) %>%
    mutate(strata = cut(y, breaks = n.bin)) %>%
    group_by(strata, month) %>%
    summarise(bf.mean = mean(BF)) 
  
  # rolling average for fill colour
  bf.df2 <- bf.df %>% 
    group_by(strata) %>%
    arrange(month) %>%
    mutate(Roll.mean = roll_mean(bf.mean, n=2, fill=0))
  
  
  #Plot
  bf.ridge <- ggplot(data= bf.df2, 
                     aes(x= month,y= strata,height = bf.mean, group = strata, fill = Roll.mean)) +
    geom_density_ridges_gradient(color="#0000000F", stat = "identity", scale=scale) +
    scale_fill_gradientn(colors = c.string,
                         limits = c(0,1),
                         name = "Mean\nProbability\nBirth") +
    scale_x_continuous(breaks = 1:12, labels=substring(month.abb, 1, 1),
                       expand = c(0,0))+
    scale_y_discrete(expand=c(0,0)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour="grey96")
    )
  
  return(bf.ridge )
}

#### Data ####
ptr.dbl.imp <- resRasterLoad("ptr", "DBL_3",T,  mod.out.dir)
##Fix names 
c.names <- names(ptr.dbl.imp)
fx.names <- gsub("ptr", "afb", c.names)
names(ptr.dbl.imp) <- fx.names
mic.dbl.imp <- resRasterLoad("mic", "DBL_3",T,  mod.out.dir)
mol.dbl.imp <- resRasterLoad("mol", "DBL_3",T,  mod.out.dir)
## adding in the dead layer for MIC
names(mic.dbl.imp)
mic.blnk <- calc(mic.dbl.imp[[1]], function(x){ifelse(is.na(x), NA, 0)})
mic.empty <- stack(mic.blnk, mic.blnk, mic.blnk)
names(mic.empty) <- c("mic1.mic2", "mic9.mic.10", "mic7.mic8")
mic.stk <-stack(mic.empty, mic.dbl.imp)
mic.fix <- mic.stk[[mixedorder(names(mic.stk))]]

ptr.I <- list(calc(ptr.dbl.imp, mean), ptr.dbl.imp)
mol.I <- list(calc(mol.dbl.imp, mean), mol.dbl.imp)
mic.I <- list(calc(mic.fix, mean), mic.fix)

#### plots ####
c1 <- colorRampPalette(c("#CFD8DC", "green4"))
c2 <- colorRampPalette(c("#CFD8DC", "dodgerblue2"))
c3 <- colorRampPalette(c("#CFD8DC", "darkorange2"))
n.ridge <- 48  # Number of ridgelines. Mess with scale below as well.
w.ridge <- .5 # Width of ridge plot compared to map

ptr.P <- BFgplotP(x = ptr.I, c.string = c1(5))
mol.P <- BFgplotP(x = mol.I, c.string = c2(5))
mic.P <- BFgplotP(x = mic.I, c.string = c3(5))

n.ridge <- 48  # Number of ridgelines. Mess with scale below as well.
w.ridge <- .5 # Width of ridge plot compared to map

rP.ptr <- BFridgeP(x = ptr.I,
                 n.bin = n.ridge,
                 scale = 2,
                 crop.extent = Africa.ext,
                 c.string = c1(5))
rP.mol <- BFridgeP(mol.I,
                 n.bin = n.ridge,scale = 2,
                 crop.extent = Africa.ext, 
                 c.string = c2(5))
rP.mic <- BFridgeP(mic.I,
                 n.bin = n.ridge,scale = 2,
                 crop.extent = Africa.ext, 
                 c.string = c3(5))
#### Post Processing ####
map.list <- list(ptr.P, mol.P, mic.P)
ridge.list <- list(rP.ptr, rP.mol, rP.mic)

map.grob <- lapply(map.list, function(x){ggplotGrob(x + guides(fill='none') +
                                                      theme(plot.margin=unit(rep(0,4), "cm")))})
ridge.grob <- lapply(ridge.list, function(x){ggplotGrob(x + guides(fill='none') +
                                                          theme(plot.margin=unit(rep(0,4), "cm"))+
                                                          coord_fixed(ratio = aspect_match))})

# find the height and left axis width of the risk grob

risk.margin.reset <- function(x){
  index <- x$layout$t[x$layout$name == 'panel']
  plot_height <- x$heights[index]
  index <- x$layout$l[x$layout$name == 'axis-l']
  axis_width <- x$widths[index]
  index <- x$layout$l[x$layout$name == 'axis-r']
  x$widths[index] <- unit(0.5, 'cm')
  return(x)
}

risk.list <- lapply(map.grob, risk.margin.reset)
# set the height the same, the left axis width, and the right
# axis to the same as the left axis for the map, scaled so that
# the widths are maintained (when laying out with grid.arrange)
# the widths are applied to the full objects, so we want everything
# to be a scaled version of the other (i.e. axis placement/widths)
# etc scaled perfectly

# Grab the correct axis width/plot heights from the first map panel
x <- map.grob[[1]]
index <- x$layout$t[x$layout$name == 'panel']
plot_height <- x$heights[index]
index <- x$layout$l[x$layout$name == 'axis-l']
axis_width <- x$widths[index]
index <- x$layout$l[x$layout$name == 'axis-r']
x$widths[index] <- unit(0.5, 'cm')

ridge.margit.reset <- function(x, plot_height, axis_width) {
  index <- x$layout$t[x$layout$name == 'panel']
  x$heights[index] <- unit(as.numeric(plot_height)/w.ridge, 'null')
  index <- x$layout$l[x$layout$name == 'axis-l']
  right_axis_width <- unit(
    grid::convertWidth(axis_width, unitTo = 'cm', valueOnly = TRUE) * w.ridge,
    "cm")
  x$widths[index] <- right_axis_width
  index <- x$layout$l[x$layout$name == 'axis-r']
  x$widths[index] <- unit(0.5 * w.ridge, 'cm')
  # x$respect <- TRUE
  return(x)
}

ridge.list <- lapply(ridge.grob, ridge.margit.reset, plot_height=plot_height, axis_width=axis_width)

master.list <- list(risk.list[[1]], ridge.list[[1]],
                    risk.list[[2]], ridge.list[[2]],
                    risk.list[[3]], ridge.list[[3]])
#### Argange and write out ####
ProbFig <- grid.arrange(grobs = master.list,
                   ncol = 2, widths=c(1,w.ridge))

ggsave(file.path(fig.pub,"probBirthComplete.pdf"),
       ProbFig,
       device = cairo_pdf,
       width=6,
       height=9.6,
       units = "in",
       dpi = 300)
