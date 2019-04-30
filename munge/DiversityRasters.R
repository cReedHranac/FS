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
  
#### ICUN Database####
## original used
mam <- readOGR(dsn = file.path(data.source, "MAMMTERR"),
            layer = "Mammals_Terrestrial")
all.genera <- sapply(mam$BINOMIAL,
                       function(x) strsplit(as.character(x), " ")[[1]][1])

## new downloaded 18.04.2019
mam.up <- readOGR(dsn = "D://TERRESTRIAL_MAMMALS",
                  layer = "TERRESTRIAL_MAMMALS")
all.equal(mam, mam.up)## Nope
## are genera conserved?
genera.up <- sapply(mam.up$binomial,
                    function(x) strsplit(as.character(x), " ")[[1]][1])
all.equal(all.genera, genera.up) ## NO


# checking which are fucked
mega.not.found <- mega.genera[which(mega.genera %!in% all.genera)] 
mega.not.found ## only occures in philipenes, not of interest for now
  ## Check with updated 
  afb.not.found <- mega.genera[which(mega.genera %!in% genera.up)] 
  #"Desmalopex"

micro.not.found <- micro.genera[which(micro.genera %!in% all.genera)] 
micro.not.found
  #"Paratriaenops" Monly found on mada and seychelles
  ## Check with updated 
  mic.not.found <- micro.genera[which(micro.genera %!in% genera.up)]
  #"Cistugo"   South Africa 
  #"Niumbaha"  african, extremely rare
  #"Hsunycteris" South american

molo.not.found <- molo.genera[which(molo.genera %!in% all.genera)]
molo.not.found
  ## Check with updated 
  mol.not.found <- molo.genera[which(molo.genera %!in% genera.up)]
  ##None!!
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

#### Create lists of all genera in Africa
genNamer <- function(x, key,  y= mam, rast = rf){
  ## function to create lists for all the genrea incloved
  all.genera <- sapply(mam$BINOMIAL,
                       function(x) strsplit(as.character(x), " ")[[1]][1])
  
  x.in <- y[which(all.genera %in% x),] #create index
  x.c <- crop(x.in, rast) #crop for speed
  x.names <- unique(x.c@data$BINOMIAL)
  Genus = sapply(x.names,
                 function(x) strsplit(as.character(x), " ")[[1]][1])
  Species = sapply(x.names,
                   function(x) strsplit(as.character(x), " ")[[1]][2])
  Class = rep(key, length(x.names))
  names.df <- as.data.frame(cbind(Genus, Species, Class))
  
  return(names.df)
  
}
genNamer.up <- function(x, key,  y= mam.up, rast = rf){
  ## function to create lists for all the genrea incloved
  ## new version since format has changed
  
  ## subset
  x.in <- y[which(y$genus %in% x),] #create index
  
  ## Cropping 
  x.c <- crop(x.in, rasterToPolygons(rast)) #crop for speed
  names.unq <- unique(x.c$binomial)
  Genus = sapply(names.unq,
                 function(x) strsplit(as.character(x), " ")[[1]][1])
  Species = sapply(names.unq,
                   function(x) strsplit(as.character(x), " ")[[1]][2])
  Class = rep(key, length(names.unq))
  names.df <- as.data.frame(cbind(Genus, Species, Class))
  
  return(names.df)
  
}
betaNator.up <- function(x,y = mam.up, rast= rf){
  ##Function for creating beta-diverstity rasters given a list of 
  ## genus names and a raster to fit it to. 
  # subset
  x.in <- y[which(y$genus %in% x),] #create index
  
  ## Cropping 
  x.c <- crop(x.in, rasterToPolygons(rast)) #crop for speed
  x.out <- rasterize(x.c, rast ,field = "binomial", fun='count', background=0)
  return(x.out)
}

afb.names <- genNamer(x= mega.genera, key = "afb")
mol.names <- genNamer(x = molo.genera, key = "mol")
mic.names <- genNamer(x = micro.genera, key = "mic")

AfricanBats <- bind_rows(afb.names, mol.names, mic.names)
write.csv(x = AfricanBats, 
          file = "data/AfricanBatSpecies.csv", 
          row.names = F)
##Update with new lists 
afb.names.up <- genNamer.up(x= mega.genera, y= mam.up, key = "afb")
mic.names.up <- genNamer.up(x= micro.genera, y= mam.up, key = "mic")
mol.names.up <- genNamer.up(x= molo.genera, y= mam.up, key = "mol")

AB.out <- bind_rows(afb.names.up,
                    mic.names.up,
                    mol.names.up)
write.csv(x = AB.out, 
          file = "data/AfricanBatSpeciesUPDATE.csv", 
          row.names = F)

afb.dist <- betaNator.up(mega.genera)
mic.dist <- betaNator.up(micro.genera)
mol.dist <- betaNator.up(molo.genera)


#### How different are these from the original ones I was using?
div.stk <- do.call(stack,
                   lapply(file.path(clean.dir, list.files(clean.dir, pattern = "*.div.tif")[c(2,3,5)]),
                          raster))
plot(afb.dist- div.stk$ptr.div)
plot(mic.dist - div.stk$mic.div)
plot(mol.dist - div.stk$mol.div)

writeRaster(afb.dist, 
            file.path(clean.dir, "ptr.div.tif"),
                      format = "GTiff",
                      overwrite = T )
writeRaster(mic.dist, 
            file.path(clean.dir, "mic.div.tif"),
            format = "GTiff",
            overwrite = T )
writeRaster(mol.dist, 
            file.path(clean.dir, "mol.div.tif"),
            format = "GTiff",
            overwrite = T )


#### Total mammalian diversity list ####
#all mammals 
all.mam <- sapply(mam$BINOMIAL,
                     function(x) strsplit(as.character(x), " ")[[1]][1])
mam.beta <- betaNator(all.mam)

### Clean ####
c.list <- c(mega, mega.genera, molo, molo.genera, micro, micro.genera, chiro, 
            chiro.genera, mam, mega.not.found, micro.not.found, molo.not.found)
rm(c.list)
