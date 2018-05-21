###############################
#### Figure 3              ####
###############################
source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)

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
BFgplot <- function(x, afr = afr.poly){
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
    scale_fill_gradientn(colors = c("#CFD8DC","#BDBDBD", "#00BCD4", "#0097A7"), 
                        limits = c(0,max(sum.df$Number)),
                        name = "Mean \nNumber \nBirthing")+

    #create african continent background
    geom_polygon(data = fortify(afr.poly),
                 aes(long, lat, group = group), 
                 colour = "#212121",
                 fill = NA) +
    aes(x=long, y=lat) +
    
    # limit coordinates
    coord_fixed(xlim = Africa.ext[1:2], ylim = Africa.ext[3:4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
    
    #Extras
    theme_bw() + 
    theme( axis.title = element_blank()) +
    scale_y_continuous(expand = c(0,0))
  
  return(g.plot)
}
BFridge <- function(x, n.bin, crop.extent = Africa.ext, scale){
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
    scale_fill_gradientn(colors = c("#CFD8DC","#BDBDBD", "#00BCD4", "#0097A7"),
                        limits = c(0,max(bf.df2$bf.mean)),
                        name = "Mean \nBirth \nForce") +
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
  
  return(bf.ridge )
}

#### Left Pannels ####
afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")

ptr.sum <- sumGen(model.string = "ptr.dbl.imp")
mol.sum <- sumGen("mol.dbl.imp")
mic.sum <- sumGen("mic.dbl.imp")


ptr.BF <- BFgplot(x = ptr.sum)
mol.BF <- BFgplot(x = mol.sum)
mic.BF <- BFgplot(x = mic.sum)

# ggsave("figures/fig3_A.png",
#        ptr.BF,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in")
# ggsave("figures/fig3_B.png",
#        mol.BF,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in")
# ggsave("figures/fig3_C.png",
#        mic.BF,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in")


#### Pannel 2 (Right) ####
n.ridge <- 48  # Number of ridgelines. Mess with scale below as well.
w.ridge <- .5 # Width of ridge plot compared to map

r.ptr <- BFridge(x = ptr.sum, n.bin = n.ridge,scale = 2,  crop.extent = Africa.ext)
r.mic <- BFridge(mic.sum, n.bin = n.ridge,scale = 2,  crop.extent = Africa.ext)
r.mol <- BFridge(mol.sum, n.bin = n.ridge,scale = 2,  crop.extent = Africa.ext)

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

png("figures/fig3_big.png", width=600, height=750, units = "mm", res = 300)
grid.arrange(grobs = master.list,
             ncol = 2, widths=c(1,w.ridge))
dev.off()
