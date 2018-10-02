##############################################################
#### Figure 1 Script                                      #### 
##############################################################

source("R/helperFunctions.R")

library(ggplot2); library(dplyr); library(data.table); library(gridExtra)
# devtools::install_github("sckott/rphylopic")
library(rphylopic);library(rgdal);library(rgeos);library(raster)
library(zoo);library(grid)

# Africa extent to use. This is reasonably tight around the data
Africa.ext <- c(-18, 47, -36, 16)

#### Figure 1 ####
#3 pannel figure, 1 pannel, map w/ outbreak locations (symbols for an/hum), 
# and 2 outbreak time lines, 1 entire history, one collapsed year
#### data ####
an.ob <- fread("data/AnOutbreakDB.csv")
an.sub <- an.ob %>%
  dplyr::select(Outbreak_ID, Org.smp, Year.Start, Month.Start, x, y) 
an.sub$Org.smp[which(an.sub$Org.smp=="Mops condylurus and Chaerephon pumilus")] <- "bat"

# an.ob <- fread(file.path(clean.dir, "annOB.PPM.csv")) #animal points
# an.an <- fread(file.path(data.source, "vIndex", "Animal_Index_Jan2018.csv")) #notes on outbreaks
# an.sub <- an.an %>%
#   dplyr::select(OUTBREAK_ID, Org.smp, Year.Start)
# names(an.sub)[1] <- names(an.ob)[1]
# an.full <- an.sub %>%
#   inner_join(an.ob,an.sub, by = "Outbreak_ID")

hum.ob <- fread("data/HumOutbreakDB.csv")
hum.sub <- hum.ob %>%
  dplyr::select(Outbreak_ID, Year.Start, Month.Start, x, y) %>%
  mutate( Org.smp = "human")

# hum.ob <- fread(file.path(clean.dir, "humOB.PPM.csv"))
# hum.an <- fread(file.path(data.source, "vIndex", "Human_Index_12_2_18.csv"))
# hum.sub <- hum.an %>%
#   dplyr::select(Outbreak_ID, Year.Start)
# hum.full <- hum.sub %>%
#   inner_join(hum.ob, hum.sub, by = "Outbreak_ID")
# 
# 
# hum.full$Org.smp <- "human"
# hum.s <- hum.full %>%
#   dplyr::select(which(names(hum.full) %in% names(an.full)))

ob.full <- bind_rows(hum.sub, an.sub)
colnames(ob.full)[which(colnames(ob.full) %in% c("x", "y"))] <- c( "long", "lat")
ob.full$Org.smp <- as.factor(ob.full$Org.smp)

z <- list()
for(i in 1:nrow(ob.full)){
  ifelse(ob.full[i,"Org.smp"] == "human", z[[i]] <- "human",z[[i]] <-  "animal")  
}

## name id ##
bat <- 	104257
chimpanzee <- 109082
duiker <-  2478895
gorilla <-  763767
human <- 	109086

critters <- c(bat, chimpanzee, duiker, gorilla, human)
critters.id <- lapply(critters, ubio_get)
critter.names <- lapply(critters.id, name_images) ## Needs manual attention
c.names <- c(critter.names[[1]]$same[[1]]$uid, 
             critter.names[[2]]$same[[1]]$uid,
             critter.names[[3]]$supertaxa[[1]]$uid,
             critter.names[[4]]$same[[1]]$uid,
             critter.names[[5]]$same[[1]]$uid)
critter.img <- lapply(c.names, image_data, size = 64)
c.img <- list()
for(i in 1:length(c.names)){
  c.img[[i]] <- critter.img[[i]][[1]]
}

## create vector of images
img.list <- list()
for(i in 1:nrow(ob.full)){
  j <- match(ob.full[i,"Org.smp"], levels(ob.full$Org.smp)) 
  img.list[[i]] <- c.img[[j]]
}

## create color vector
cz <- c("darkorange2", "black", "green4", "dodgerblue2", "red")
col.list <- list()
for(i in 1:nrow(ob.full)){
  j <- match(ob.full[i,"Org.smp"], levels(ob.full$Org.smp)) 
  col.list[[i]] <- cz[[j]]
}

## lets create a bounding box polygon
#bounds
bbx <- c(left = 12.5,
         bottom= -.25,
         right= 15.5,
         top= 1.5)

x <- c(bbx["left"], bbx["left"], bbx["right"], bbx["right"])
y <- c(bbx["bottom"], bbx["top"], bbx["top"], bbx["bottom"])
df <- data.frame(x, y)

#### Pannel 1 Maps ####
## accessory layers
region_colour <- '#CFD8DC'

afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")
rf.poly <- rasterToPolygons(raster(file.path(data.source, "cropMask.tif")),
                            fun = function(x){x==1}, dissolve = T)
bkg <- theme_bw() + theme( 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank()) 

ob.plot <- ggplot() +
  geom_polygon(data = fortify(rf.poly),
               aes(long, lat, group = group),
               colour = "#212121",
               fill = region_colour,
               size = 0.1)+
  geom_polygon(data = fortify(afr.poly),
               aes(long, lat, group = group), 
               colour = "#212121",
               fill = "NA") +
  coord_fixed(xlim = Africa.ext[1:2], ylim = Africa.ext[3:4]) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-20, 50, by=10)) +
  theme_bw()


for( i in 1:nrow(ob.full)){
  if(ob.full$Org.smp[i] == "human"){
    j <- 3
  } else{ifelse(ob.full$Org.smp[i] == "bat",  j <- 1.9,  j <- 2.5)}
  ob.plot <- ob.plot + 
    add_phylopic(img.list[[i]], .5, 
                 ob.full$long[i],
                 ob.full$lat[i],
                 ysize = j,
                 color = col.list[[i]]) 
  
}

ob.plot <- ob.plot +
  geom_polygon(aes(x=x, y=y), data=df, color = "red", alpha=.5)
  
## good enough for now...


### Insert plot
insert.ext <- c(12.5, 15.5,1.5,-.25)
ob.insert <- ggplot()+
  geom_polygon(data = fortify(afr.poly),
               aes(long, lat, group = group), 
               colour = "#212121",
               fill = region_colour) +
  coord_fixed(xlim = insert.ext[1:2],ylim = insert.ext[3:4]) +
  theme_bw() + theme(axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title = element_blank())

for( i in 1:nrow(ob.full)){
  if(ob.full$Org.smp[i] == "human"){
    j <- .3
  } else{ifelse(ob.full$Org.smp[i] == "bat",  j <- .2,  j <- .25)}
  ob.insert <- ob.insert + 
    add_phylopic(img.list[[i]], .5, 
                 ob.full$long[i],
                 ob.full$lat[i],
                 ysize = j, ##humans need bigger, bats smaller
                 color = col.list[[i]]) 
}
bkg.insert <- theme(
  plot.margin = unit(c(0,0,0,0,0), "mm"),
  panel.border = element_rect(colour = "red", fill=NA, size=1))

ob.insert +bkg.insert ## That's pretty good... 


map.with.insert <- ob.plot + bkg +
  annotation_custom(grob = ggplotGrob(ob.insert+bkg.insert),
                    xmin = -17,
                    xmax = 17.5, 
                    ymin = -36,
                    ymax = -14.5)

# ggsave("figures/fig1_A.png",
#        map.with.insert,
#        device = "png",
#        width = 210,
#        height = 165,
#        units = "mm", 
#        dpi = 300)

#### Pannel 2 Time line ####
ob.T <- ob.full
ob.T$Date <- as.yearmon(paste(ob.T$Year.Start, ob.T$Month.Start), "%Y %m")
ob.a <- ob.T %>%
  dplyr::arrange(Date)
## Functions to sort the phylo pics accordingly addapted from https://stackoverflow.com/questions/27637455/display-custom-image-as-geom-point
RlogoGrob <- function(x, y, size, img, col, alpha) {
  mat = rphylopic:::recolor_phylopic(img, alpha, col)
  aspratio <- ncol(mat)/nrow(mat)
  rasterGrob(x = x, y = y, image = mat, default.units = "native", 
             height = size*2.5)#, 
             #width = size*aspratio*2.5)
}

GeomRlogo <- ggproto("GeomRlogo", Geom, draw_panel = function(data, panel_scales, 
                                                              coord, img, colour, alpha, na.rm = FALSE) {
  coords <- coord$transform(data, panel_scales)
  ggplot2:::ggname("geom_Rlogo", RlogoGrob(coords$x, coords$y, coords$size, 
                                           img, colour, alpha))
}, non_missing_aes = c("Rlogo", "size"), required_aes = c("x", "y"), default_aes = aes(size = 0.05), 
icon = function(.) {
}, desc_params = list(), seealso = list(geom_point = GeomPoint$desc), 
examples = function(.) {
})

geom_Rlogo <- function(mapping = NULL, data = NULL, stat = "identity", 
                       position = "identity", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, 
                       img,size,  ...) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomRlogo, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(na.rm = na.rm, img = img, size = size,  ...))
}

## Plot 
y_vals <- 5:1 #seq(20,20*5,by=20)
s <- 0.05
ob.plot <- ob.a %>% mutate(y = y_vals[as.numeric(as.factor(Org.smp))])

g.time <- ggplot(data = ob.plot, aes(x= Date, y = y)) +
  geom_point(data = expand.grid(Date=filter(ob.a, Org.smp == "human") %>% dplyr::pull(Date), y=y_vals), ## Grey points
             aes(x = Date, y= y, size = 1, alpha = .5),
             color = "grey70", shape=20,
             show.legend = F)+
  geom_point(aes(x = Date, y= y, color = Org.smp, alpha = .5), size = 4,
             show.legend = F)+ 
  geom_segment(data=data.frame(y=1:5), aes(x = 1975, y = y_vals, xend = 2018.5, yend = y_vals),
               arrow = arrow(length =  unit(x = 0.2,units = 'cm'),type = 'closed')) +
  scale_x_yearmon(format = "%Y", n = 10, limits = c(1973.5, 2019), expand=c(0,0)) +
  scale_y_continuous(limits=c(0.5,5.5), expand=c(0,0)) +
  geom_Rlogo(aes(x, y), img=c.img[[1]], alpha=1, col="darkorange2", size = s, data=data.frame(x=1975, y=y_vals[1], Org.smp = "bat" )) +
  geom_Rlogo(aes(x, y), img=c.img[[2]], alpha=1, col="black", size = s, data=data.frame(x=1975, y=y_vals[2], Org.smp = "chimpanzee" )) +
  geom_Rlogo(aes(x, y), img=c.img[[3]], alpha=1, col='green4', size = s*0.15/0.125, data=data.frame(x=1975, y=y_vals[3], Org.smp = "duiker" )) +
  geom_Rlogo(aes(x, y), img=c.img[[4]], alpha=1, col='dodgerblue2', size = s*0.15/0.125, data=data.frame(x=1975, y=y_vals[4], Org.smp = "gorilla" )) +
  geom_Rlogo(aes(x, y), img=c.img[[5]], alpha=1, col='red', size = s*0.18/0.125, data=data.frame(x=1975, y=y_vals[5], Org.smp = "human" )) +
  scale_colour_manual(values = c(bat = 'darkorange2',
                                 chimpanzee ='black',
                                 duiker = 'green4',
                                 gorilla = 'dodgerblue2',
                                 human = 'red')) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor.y = element_blank())

g.time

# ggsave("figures/fig1_B.png",
#        g.time,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in",
#        dpi = 300)

#### Pannel 3 Bar with Smooth ####
### If we can start this at setptember we can likely get the bimodal distribution
### to show up better. Do two different smooths for human/nonhumans and maybe one for total?
ob.hist <- ob.a

org=data.frame(Org.smp=levels(as.factor(ob.hist$Org.smp)),
               Org=c("Chiroptera",
                     "Non-Volant Mammals",
                     "Non-Volant Mammals",
                     "Non-Volant Mammals",
                     "Human"))

ob.hist <- inner_join(ob.hist, org, "Org.smp")

ob.hist$Month <- factor(format(ob.hist$Date, "%b"), 
                        levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
ob.hist$Org <- factor(ob.hist$Org, levels(ob.hist$Org)[c(1,3,2)])
g.bar <- ggplot(data=ob.hist,aes(x=Month, fill=Org.smp))+
  geom_bar(show.legend = F) +
  scale_fill_manual(values = c(bat = 'darkorange2',
                                 chimpanzee ='black',
                                 duiker = 'green4',
                                 gorilla = 'dodgerblue2',
                                 human = 'red')) +
  scale_x_discrete(labels=substring(month.abb, 1, 1)) +
  facet_wrap(~Org, nrow = 3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,7.2)) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  )

g.bar


# ggsave("figures/fig1_C.png",
#        g.bar,
#        device = "png",
#        width = 5,
#        height = 5,
#        units = "in",
#        dpi = 300)

### put them together

fig1.complete <- grid.arrange(map.with.insert, g.time, g.bar,
             widths = c(2.5,1), heights = c(2.5,1),
             layout_matrix = rbind(c(1,3),
                                   c(2,3)))

ggsave("figures/Fig1Complete.eps",
       fig1.complete,
       device = cairo_ps, 
       width = 7.5,
       height = 6,
       units = "in",
       dpi = 300)
