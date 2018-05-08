##########################
#### Figure 4         ####
##########################

# feel free to specify the path to FS (including the FS folder) in base.path.
# if not set, it'll use what Reed uses (D:\Dropbox\FS or ~/Dropbox/FS)
#base.path <- "~/data/Dropbox/FS"

source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

# Africa extent to use. This is reasonably tight around the data
Africa.ext <- c(-18, 47, -36, 16)

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
    theme_bw() +
    theme(axis.title = element_blank())

  return(g.plot)
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

#### Alternative human models for ridges ####

hum.NoAn <- spatHandler("humNoAn") # Go with this one for ridges

n.ridge <- 96  # Number of ridgelines. Mess with scale below as well.
w.ridge <- 0.5 # Width of ridge plot compared to map

afr.ridge <- ERridge(hum.NoAn, n.bin = n.ridge, scale = 2, crop.extent = Africa.ext )

# Arrange the map and ridge on the same plot

# first work out aspect ratio of map
width_height <- diff(Africa.ext)[c(1,3)]
aspect_map <- width_height[1] / width_height[2]

# now aspect ratio of ridge plot (this assumes top bin isn't too large)
aspect_ridge <- 12 / n.ridge

# and fixup the aspect ratio of ridge to that of map
aspect_match <- aspect_ridge / aspect_map / w.ridge

risk.grob <- ggplotGrob(risk.plot + guides(fill='none') +
                          theme(plot.margin=unit(rep(0,4), "cm")))
ridge.grob <- ggplotGrob(afr.ridge + guides(fill='none') +
                           theme(plot.margin=unit(rep(0,4), "cm")) +
                           coord_fixed(ratio = aspect_match)) # Alternatively see the respect=TRUE below

# find the height and left axis width of the risk grob
index <- risk.grob$layout$t[risk.grob$layout$name == 'panel']
plot_height <- risk.grob$heights[index]
index <- risk.grob$layout$l[risk.grob$layout$name == 'axis-l']
axis_width <- risk.grob$widths[index]
index <- risk.grob$layout$l[risk.grob$layout$name == 'axis-r']
risk.grob$widths[index] <- unit(0.5, 'cm')

# set the height the same, the left axis width, and the right
# axis to the same as the left axis for the map, scaled so that
# the widths are maintained (when laying out with grid.arrange)
# the widths are applied to the full objects, so we want everything
# to be a scaled version of the other (i.e. axis placement/widths)
# etc scaled perfectly
index <- ridge.grob$layout$t[ridge.grob$layout$name == 'panel']
ridge.grob$heights[index] <- unit(as.numeric(plot_height)/w.ridge, 'null')
index <- ridge.grob$layout$l[ridge.grob$layout$name == 'axis-l']
right_axis_width <- unit(
  grid::convertWidth(axis_width, unitTo = 'cm', valueOnly = TRUE) * w.ridge,
  "cm")
ridge.grob$widths[index] <- right_axis_width
index <- ridge.grob$layout$l[ridge.grob$layout$name == 'axis-r']
ridge.grob$widths[index] <- unit(0.5 * w.ridge, 'cm')
#ridge.grob$respect <- TRUE # This is equivalent to coord_fixed() - i.e. it sets a constant aspect ratio.

# arrange them side by side
png("figures/fig4_Complete.png", width=800, height=430)
grid.arrange(risk.grob,
             ridge.grob,
             ncol=2, widths=c(1,w.ridge))
dev.off()



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

central.ridge <- ERridge(hum.NoAn, n.bin = 30, central.africa)
western.ridge <- ERridge(hum.NoAn, n.bin = 30, west.africa)

#### Animal Risk Plot ####
ann.mean <- spatHandler("ann")
ann.risk.plot <- ERgplot(ann.mean)
ann.ridge <- ERridge(ann.mean, n.bin = n.ridge, scale = 2, crop.extent = Africa.ext )
# and fixup the aspect ratio of ridge to that of map
aspect_match <- aspect_ridge / aspect_map / w.ridge

risk.grob <- ggplotGrob(ann.risk.plot + guides(fill='none') +
                          theme(plot.margin=unit(rep(0,4), "cm")))
ridge.grob <- ggplotGrob(ann.ridge + guides(fill='none') +
                           theme(plot.margin=unit(rep(0,4), "cm")) +
                           coord_fixed(ratio = aspect_match))

# find the height and left axis width of the risk grob
index <- risk.grob$layout$t[risk.grob$layout$name == 'panel']
plot_height <- risk.grob$heights[index]
index <- risk.grob$layout$l[risk.grob$layout$name == 'axis-l']
axis_width <- risk.grob$widths[index]
index <- risk.grob$layout$l[risk.grob$layout$name == 'axis-r']
risk.grob$widths[index] <- unit(0.5, 'cm')

# set the height the same, the left axis width, and the right
# axis to the same as the left axis for the map, scaled so that
# the widths are maintained (when laying out with grid.arrange)
# the widths are applied to the full objects, so we want everything
# to be a scaled version of the other (i.e. axis placement/widths)
# etc scaled perfectly
index <- ridge.grob$layout$t[ridge.grob$layout$name == 'panel']
ridge.grob$heights[index] <- unit(as.numeric(plot_height)/w.ridge, 'null')
index <- ridge.grob$layout$l[ridge.grob$layout$name == 'axis-l']
right_axis_width <- unit(
  grid::convertWidth(axis_width, unitTo = 'cm', valueOnly = TRUE) * w.ridge,
  "cm")
ridge.grob$widths[index] <- right_axis_width
index <- ridge.grob$layout$l[ridge.grob$layout$name == 'axis-r']
ridge.grob$widths[index] <- unit(0.5 * w.ridge, 'cm')

# arrange them side by side
png("figures/AnnRisk.png", width=800, height=430)
grid.arrange(risk.grob,
             ridge.grob,
             ncol=2, widths=c(1,w.ridge))
dev.off()

