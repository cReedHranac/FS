##### Helper functions for FiloSpatial #####
##### Re-Script Sept 2017              #####
############################################

#### Platform specific path assignment #####
if(.Platform$"OS.type" == "windows"){
  data.source <-file.path("D:", "Dropbox", "FS", "SourceData")
  clean.dir <- file.path("D:", "Dropbox", "FS", "Processed")
  norm.dir <- file.path("D:", "Dropbox", "FS", "Normalized")
} else{
  data.source <-file.path("~", "Dropbox", "FS", "SourceData")
  clean.dir <- file.path("~", "Dropbox", "FS", "Processed")
  norm.dir <- file.path("~", "Dropbox", "FS", "Normalized")
}

#### Opperators ####
## The %!in% opperator 
'%!in%' <- function(x,y)!('%in%'(x,y))

#### Loading Funs ####
## occLoad_fun
occLoad <- function(taxon.class){
  ##Function for loading the occurence dataframes from cache 
  #Arguments:
  #taxon.class <- 3 letter sub set that corisponds with the breeding class
  # options built in are "ptr", "mol", "mic" (and "all" not currently supported)
  
  ## Check input
  if(taxon.class %!in% c("ptr", "mol", "mic")){
    cat("Taxon class is not recognized try again")
  }
  if(taxon.class == "ptr"){
    assign("ptr.occ", read.csv(file.path(clean.dir, "ptrOcc.csv")),
           envir = .GlobalEnv)
  }
  if(taxon.class == "mol"){
    assign("mol.occ", read.csv(file.path(clean.dir, "molOcc.csv")),
           envir = .GlobalEnv)
  }
  if(taxon.class == "mic"){
    assign("mic.occ", read.csv(file.path(clean.dir, "micOcc.csv")),
           envir = .GlobalEnv)
  }
}

staticStkLoad <- function(path.to.dir = clean.dir){
  ##Function the load the static covarriates from storage.
  ## Loads them to the global env with name 'static.stk'
  #Arguments:
    #path.to.dir <- generally speaking will be clean.dir but it coule be the 
    #any other if needed (and as long as the names are consistant)
  #N.B. These layers will not be normalized
  require(raster)
  
  static <- c( sub.files <- paste0("bio",c(2,4,7,15),".tif"), #BClim subset
               list.files(clean.dir, pattern = "*_sum.tif"),  #Diversity Rasters
               "LandCover.tif")                               #LandCover
  stk <- do.call(stack, lapply(file.path(clean.dir, static), raster))
  assign("static.stk", stk, envir = .GlobalEnv)
}

#### Logs and Directories ####
outDirGen <- function(run.id, mod.id){
  out.dir <- file.path("ModOut", run.id, mod.id)  
  if(!file.exists(out.dir)){
    dir.create(file.path(outdir),recursive=T, showWarnings = F)
  } 
  return(out.dir)
}

log.path <- function(path.to.dir = out.dir, run.id, mod.id){
  out.log <- file.path(out.dir, paste0("TextLog_",run.id,"_",mod.id,".txt"))
}


#### Pre-Processing Functions ####
listGen <- function(occ.db){
  ### Function to generate the list of names required for an lapply
  ## Arguments:
  ##  occ.db: product of doubleMonth() or 
  occ.l <- names(occ.db[3:(length(names(occ.db))-1)])
  return(occ.l)
}

envPCA <- function(env, mod.id = NULL, path.to.dir = NULL){
  ### Function to reduce env dimentionality prior to biomod modeling
  ### Preforms PCA on all varriables aside from LandCover and adds LandCover back
  ## Arguments:
  ## env: RasterStack (generally produced from stkRevolver) to preform operation on
  ## mod.id: myRespName for out put file needs
  ## file.out: outdir for out put file needs
  require(RStoolbox)
  cont.env <- env[[-which(names(env)=="LandCover")]]
  env.pca <- rasterPCA(cont.env, nComp = 6, spca = F)
  env.stk <- stack(unstack(env.pca$map))
  env.out <- addLayer(env.stk, env[[which(names(env)=="LandCover")]])
  if(!is.null(name) || !is.null(file.out)){
    capture.output(env.pca$model$loadings, 
                   file = file.path(path.to.dir, paste0("PCA_LOAD_",mod.id,".txt")))
    capture.output(env.pca,
                   file = file.path(path.to.dir, paste0("PCA_LOAD_",name,".txt")), append = T)
  }
  return(env.out)
}


