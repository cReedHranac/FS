##############################
#### Figure 2 Script      ####
##############################
source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra)
### Data #### 
library(rphylopic)
bat <- 	104257
bat.u <- ubio_get(bat)
bat.n <- name_images(bat.u)
bat.symbol <- image_data(bat.n$same[[1]]$uid, size = 64)[[1]]

bi <- fread("data/afrBatBirthDB.csv")
bi$Class <- factor(bi$Class,
                   labels= c("Pteropodiae", "Molosidae", "Non-Molosid Microbats"))
bi$Start <- factor(bi$Start,
                   labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
bi$long <- bi$Lon
bi$lat <- bi$Lat



cols <-  c("green4", "dodgerblue2", "darkorange2")
col.list <- list()
for (i in 1:nrow(bi)) {
  j <- match(bi$Class[i], levels(bi$Class))
  col.list[[i]] <- cols[[j]]
}

afr.poly <- readOGR(dsn = file.path(data.source, "Africa"),
                    layer = "AfricanCountires")
rf.poly <- rasterToPolygons(raster(file.path(data.source, "cropMask.tif")),
                            fun = function(x){x==1}, dissolve = T)
bkg <- theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  plot.title = element_text(hjust = 0.5))



bi.t <- bi %>%
  group_by(Start) %>%
  summarise(num = n())

#### Pannel 1 Map ####
bi.plot <- ggplot()+ bkg +
  geom_polygon(data = fortify(afr.poly),
               aes(long, lat, group = group), 
               colour = "grey20",
               alpha = .25) +
  geom_polygon(data = fortify(rf.poly),
               aes(long, lat, group = group),
               colour = "white", 
               alpha = .25,
               fill = "cornsilk")+
  coord_fixed(xlim = c(-18, 49),ylim = c(-36, 15))

for(i in 1:nrow(bi)){
  bi.plot <- bi.plot +
    add_phylopic(bat.symbol, .65,
                 x = bi$long[i],
                 y = bi$lat[i],
                 ysize = 1.7, 
                 color = col.list[[i]])
}

bi.plot ## looks pretty good. still need to sort the ledgends though

ggsave("figures/fig2_A.png",
       bi.plot,
       device = "png",
       width = 5,
       height = 5,
       units = "in")

#### Pannel 2 Bar with Smooth ####
bi.bar <- ggplot(data = bi)+
  geom_bar(aes(x=Start, fill = Class),
           show.legend = F) +
  scale_fill_manual(values = c("green4", "dodgerblue2", "darkorange2"))+
  facet_wrap(~Class, ncol = 1, scales = "free_y") +
  theme_bw() + 
  theme(
    axis.ticks.x = element_blank(),
    axis.title = element_blank())

bi.bar

ggsave("figures/fig2_B.png",
       bi.plot,
       device = "png",
       width = 5,
       height = 5,
       units = "in")
