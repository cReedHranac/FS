#######################################
## Diversity Maps
#######################################
## Bat diversity by class
##Reading in the mammallian shapefiles
#Creating an empty raster for template use
# b1 <- raster("data/AfrPrecip/Precip25/AfrPrecip25_1.tif")
# as.1 <- function(x){
#   if(!is.na(x)){return(1)}else{return(NA)}
# }
# afr.blnk <- calc(b1, as.1)
# rm(b1)
#mam <- readOGR(dsn = "data/MAMMTERR",
#               layer = "Mammals_Terrestrial")
# Using taxize to get all bat genera
#  library(taxize)
#  chiro <- itis_downstream( tsns = 179985, downto = 'Genus')
#   #getting genus names
#   chiro.genera <- as.character(chiro$taxonname)
#   # writing out
#   write.csv(chiro.genera, "data/Dataframes/chiro_genera.csv", row.names = F)
# # using taxize to get a list of genera for Pteropodidae (megachiroptera)
# mega <- itis_downstream( tsns = 180029 ,downto = 'Genus')
#   #getting genus names
#   mega.genera <- as.character(mega$taxonname)
#   #wrting out
#   write.csv(mega.genera, "data/Dataframes/mega_genera.csv", row.names = F)
# # taxize for Molossidae
# molo <- itis_downstream(tsns = 180077, 'Genus')
#   #getting genus names
#   molo.genera <- as.character(molo$taxonname)
#   #writing out
#   write.csv(molo.genera, "data/Dataframes/molo.genera.csv", row.names=F)
# # taxize for Ynagochiroptera (microchiroptera)
# micro <- itis_downstream(tsns = 945841, 'Genus')
#   #getttign genera names
#   micro.genera <- as.character(micro$taxonname)
#   #writing out
#   write.csv(micro.genera, "data/Dataframes/micro_genera.csv", row.names = F)

# chiro.genera <- read.csv("data/Dataframes/chiro_genera.csv")
# chiro.genera <- chiro.genera[,1]
# mega.genera <- read.csv("data/Dataframes/mega_genera.csv")
# mega.genera <- mega.genera[,1]
# micro.genera <- read.csv("data/Dataframes/micro_genera.csv")
# micro.genera <- micro.genera[,1]
# molo.genera <- read.csv("data/Dataframes/molo.genera.csv")
# molo.genera <- molo.genera[,1]
# 
# # finding all genera names from the Mammterr files
#   chiro <- itis_downstream( tsns = 179985, downto = 'Species')
# 
# all.genera <- sapply(mam$BINOMIAL,
#                       function(x) strsplit(as.character(x), " ")[[1]][1])
# 
#   # Checking to make sure that all genera are accounted for
#   chiro.match <- which(chiro.genera %in% all.genera)
#   chiro.matched <- chiro.genera[chiro.match]
#   not.in <- setdiff(chiro.genera, chiro.matched)
#     # Addressing
#       #"Desmalopex" Phillopenese only
#       #"Paratriaenops" used to be triaenops
#       ("Triaenops" %in% all.genera) # True 
#       #"Chaerephon" distinction of Tadarida / mop
#       ("Tadarida" %in% all.genera) # True
#       # "Paremballonura"  used to be Emballonura
#       ("Emballonura" %in% all.genera) #True
#       #"Dryadonycteris" brazillian single species
#       #"Hypsugo" formerrly Pipistrellus
#       ("Pipistrellus" %in% all.genera) # True
#       #"Neoromicia" formerly either pipistrellus, Eptesicus or Vespertilio
#       ("Vespertilio" %in% all.genera) # True
#       ("Eptesicus" %in% all.genera) # True
#       #Eptesicus  formerlly Artibeus
#       ("Artibeus" %in% all.genera) # True
#       #"Parastrellus" formerlly Pipistrellus
#       #"Perimyotis" same ^
#       # "Nimbaha" formerlly glauconyteris
#       ("Glauconyteris" %in% all.genera) #False  (Very rare)
#       #"Hsunycteris" brazillian only
#       # Vampyriscus formerlly subgenera of Vampyressa
#       ("Vampyressa" %in% all.genera) #True
#     #Checking to make sure that those occurr in the list of chiro genrea 
#   check.v <- c("Triaenops", "Tadarida", "Emballonura", "Pipistrellus", "Vespertilio","Eptesicus","Artibeus", "Vampyressa")
#   (check.v %in% chiro.genera)  # all TRUE

# # find all chiropterean genera within the list
# chiro.idx <- which(all.genera %in% chiro.genera)
# 
#   mega.index <- which(all.genera %in% mega.genera)
#   micro.index <- which(all.genera %in% micro.genera)
#   molo.index <- which(all.genera %in% molo.genera)
# # get the shapefiles for each   
# chiro.shp <- mam[chiro.idx, ]
# 
#   mega.shp <- mam[mega.index,]
#   micro.shp <- mam[micro.index,]
#   molo.shp <- mam[molo.index,]
# # # remove mammterr for memory
# # rm(mam)
# 
# # get the names of species for whuch there are shapefiles
# chiro.species <- unique(chiro.shp$BINOMIAL)
#   mega.spp <- unique(mega.shp$BINOMIAL)
#   micro.spp <- unique(micro.shp$BINOMIAL)
#   molo.spp <- unique(molo.shp$BINOMIAL)
#   
#   
# rast.dist <- function(sp, mask.raster=afr.blnk){
#     #function to rasterize the species distributions
#   
#   shp <- mam[mam$BINOMIAL == sp,]
#   
#   r <- rasterize(shp, mask.raster, field=1, fun= 'first', background=0, mask=T)
#   
#   r@data@names <- shp@data$BINOMIAL
#   
#   if(r@data@min==1){return(r)}
# }
# 
# 
# #initalizing stack
# mega <- list()
# mega <- lapply(mega.spp, rast.dist)
# mega.stk <- do.call(stack, mega)
# mega.sum <- stackApply(mega.stk,1, sum)
# 
# micro <- ls()
# micro <- lapply(micro.spp, rast.dist)
#   #removing null elements
# micro <- micro[!sapply(micro, is.null)]
# micro.stk <- do.call(stack,micro)
# micro.sum <- stackApply(micro.stk,1,sum)
# 
# molo <- stack()
# molo <- lapply(molo.spp, rast.dist)
# molo.stk <- do.call(stack, molo)
# molo.sum <- stackApply(molo.stk,1,sum)
# 
# # All mammallian species across all Africa
# mam.gen <- unique(mam@data$BINOMIAL)
# all.mam <- list()
# all.mam <- lapply(mam.gen, rast.dist, stk=all.mam, stk.name = "mam")
# all.mam <- all.mam[!sapply(all.mam, is.null)]
# mam.stk <- do.call(stack, all.mam)
# mam.sum <- stackApply(mam.stk,1,sum)
# 
# writeRaster(mega.sum, "data/DiversityRasters/Mega_sum.tif")
# writeRaster(micro.sum, "data/DiversityRasters/Micro_sum.tif")
# writeRaster(molo.sum, "data/DiversityRasters/Molo_sum.tif")
# writeRaster(mam.sum, "data/DiversityRasters/Mam_sum.tif")

files <- list.files("data/DiversityRasters/")
div.rast <- lapply(paste0("data/DiversityRasters/",files), raster)
div.rast
rm(files)

## Ammendment: Create non-bat diversity layer 
library(raster)
div.stk <- stack(file.path(clean.dir, 
                           list.files(file.path(clean.dir), pattern = "*_sum")))
nb.div <- div.stk$Mam_sum - (div.stk$Mega_sum + div.stk$Micro_sum + div.stk$Molo_sum) 
writeRaster(nb.div, file.path(clean.dir, "NB_sum.tif"), format = "GTiff")
