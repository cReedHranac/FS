#######################################
## Diversity Maps
#######################################
## Bat diversity by class
##Reading in the mammallian shapefiles
#Creating an empty raster 
source("R/helperFunctions.R")

rf <- raster(file.path(data.source,"cropMask.tif"))

##Using taxize to get all bat genera
 library(taxize)
 # using taxize to get a list of genera for Pteropodidae (megachiroptera)
mega <- itis_downstream( tsns = 180029 ,downto = 'Genus')
  #getting genus names
mega.genera <- as.character(mega$taxonname)
 # taxize for Molossidae
molo <- itis_downstream(tsns = 180077, 'Genus')
  #getting genus names
molo.genera <- as.character(molo$taxonname)
# taxize for Ynagochiroptera (microchiroptera)
micro <- itis_downstream(tsns = 945841, 'Genus')
  #getttign genera names
micro.genera <- as.character(micro$taxonname)

## Double checking
chiro <- itis_downstream( tsns = 179985, downto = 'Species')
all.chiro <- as.character(chiro$taxonname)

all(mega.genera %in% all.chiro)
all(micro.genera %in% all.chiro)  
all(molo.genera %in% all.chiro)
  
### ICUN Database
mam <- readOGR(dsn = file.path(data.source, "MAMMTERR"),
            layer = "Mammals_Terrestrial")
all.genera <- sapply(mam.o$BINOMIAL,
                       function(x) strsplit(as.character(x), " ")[[1]][1])

# checking which are fucked
mega.not.found <- mega.genera[which(mega.genera %!in% all.genera)] 
mega.not.found ## only occures in philipenes, not of interest for now

micro.not.found <- micro.genera[which(micro.genera %!in% all.genera)] 
micro.not.found
  #"Paratriaenops" Monly found on mada and seychelles

molo.not.found <- molo.genera[which(molo.genera %!in% all.genera)]
molo.not.found

####NB: There are many issues found particullarly strange is the lack of Chaerephon & Mops. ####
## attempted to get the more up-to-date IUCN database and ~30,000 of the species are missing including
## a large number of the bat species.
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
#       #
#       # Vampyriscus formerlly subgenera of Vampyressa
#       ("Vampyressa" %in% all.genera) #True
#     #Checking to make sure that those occurr in the list of chiro genrea
#   check.v <- c("Triaenops", "Tadarida", "Emballonura", "Pipistrellus", "Vespertilio","Eptesicus","Artibeus", "Vampyressa")
#   (check.v %in% chiro.genera)  # all TRUE
#### Continue ####
betaNator <- function(x,y = mam, rast= rf){
  ##Function for creating beta-diverstity rasters given a list of 
  ## genus names and a raster to fit it to. 
  all.genera <- sapply(mam$BINOMIAL,
                       function(x) strsplit(as.character(x), " ")[[1]][1])
  
  x.in <- y[which(all.genera %in% x),] #create index
  x.c <- crop(x.in, rast) #crop for speed
  x.out <- rasterize(x.c, rast ,field = "BINOMIAL", fun='count', background=0)
  return(x.out)
}
ptr.beta <- betaNator(mega.genera)
mic.beta <- betaNator(micro.genera)
mol.beta <- betaNator(molo.genera)

#all mammals 
all.mam <- sapply(mam$BINOMIAL,
                     function(x) strsplit(as.character(x), " ")[[1]][1])
mam.beta <- betaNator(all.mam)

### Clean ####
c.list <- c(mega, mega.genera, molo, molo.genera, micro, micro.genera, chiro, 
            chiro.genera, mam, mega.not.found, micro.not.found, molo.not.found)
rm(list= c.list)