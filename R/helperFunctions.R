##### Helper functions for FiloSpatial #####
##### Re-Script Sept 2017              #####
############################################

#### Platform specific path assignment #####
if (!exists('base.path')) {
  if(.Platform$"OS.type" == "windows"){
    base.path = file.path("D:", "Dropbox", "FS")
  } else {
    base.path = path.expand("~/Dropbox/FS")
  }
}
data.source <-file.path(base.path, "SourceData")
clean.dir <- file.path(base.path, "Processed")
clean.dir.nov <- file.path(base.path, "novAddition")
norm.dir <- file.path(base.path, "Normalized")
mod.out.dir <- file.path(base.path, "ModOut")

#### Opperators ####
## The %!in% opperator 
'%!in%' <- function(x,y)!('%in%'(x,y))

#### Loading Funs ####

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
               list.files(clean.dir, pattern = "*.div.tif")[c(2:5)],  #Diversity Rasters
               "LandCover.tif")                               #LandCover
  stk <- do.call(stack, lapply(file.path(clean.dir, static), raster))
  assign("static.stk", stk, envir = .GlobalEnv)
}

stkRevolver <- function(month, path.to.dir = norm.dir){
  ### Function for obtaining temporal covarriates within the requested window
  ## Arguments: 
  ##  month <- month (as intiger) for the cylinder to start on
  # can also be a name item from a column header
  ##  path.to.dir <- directory to look for the layers in .
  # normalized rasters will be in "norm.dir"
  # processed but un-altered rasters will be in "clean.dir"
  # source rasters will be in "source.dir"
  
  require(raster);require(gtools)
  tempo <- c("evi", "pet", "prec", "tmean") #temporal co-varriates
  
  internalRevolver <- function(x, start.i = month){
    ##Internal function to read in the subset of rasters required, and stack  
    if("numeric" %!in% is(start.i)){
      start.i <- as.numeric(substring(start.i,4,4))
    }
    
    ##Heart of the Revolver
    v <- 1:12
    i <- ((v + start.i) - 2) %% length(v)+1
    j <- i[c(1:3,11:12)] ## number vector for 1 month
    
    ##Read in and stack
    stk <- do.call(stack, lapply(file.path(path.to.dir,
                                           mixedsort(list.files(path.to.dir,
                                                                pattern = paste0(x,"*")))[j]),
                                 raster))
    return(stk)
  }
  
  #apply internal function and stack products
  tempo.env <- do.call(stack, lapply(tempo, internalRevolver, start.i = month))
  
  return(tempo.env)
}

#### Logs and Directories ####
outDirGen <- function(mod.id, run.id){
  out.dir <- file.path("ModOut", run.id, mod.id)
  if(!file.exists(out.dir)){
    dir.create(file.path(out.dir),recursive=T, showWarnings = F)
  }
  return(out.dir)
}

log.path <- function(path.to.dir = out.dir, run.id, mod.id){
  out.log <- file.path(path.to.dir, paste0("TextLog_",run.id,"_",mod.id,".txt"))
}


#### Pre-Processing Functions ####
listGen <- function(occ.db){
  ### Function to generate the list of names required for an lapply
  ## Arguments:
  ##  occ.db: product of occLoad or occLoad2
  occ.l <- names(which(apply(occ.db[,3:ncol(occ.db)],2,sum) > 5))
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
  if(!is.null(mod.id) || !is.null(file.out)){
    capture.output(env.pca$model$loadings, 
                   file = file.path(path.to.dir, paste0("PCA_LOAD_",mod.id,".txt")))
    capture.output(env.pca,
                   file = file.path(path.to.dir, paste0("PCA_LOAD_",mod.id,".txt")), append = T)
  }
  return(env.out)
}

#### Varriable Importance Functions ####
normImp <- function(BiomodOut){
  ### Function for normalizing varriable importance scores
  ## Argument:
  ## Biomodout:  Model out to retrive scores from 
  require(dplyr)
  varImportance <- get_variables_importance(BiomodOut)
  row.n <- rownames(varImportance)
  col.n <- colnames(varImportance)
  dim(varImportance) <- dim(varImportance)[-4]   # get rid of fourth dimension
  
  list_of_VI <- lapply(1:dim(varImportance)[3], function(i,...) { varImportance[,,i] }) # Stack to list of matrix
  norm_Vi <- lapply(list_of_VI,normVar)  # Normalize variable scors
  bfdf <- as.data.frame(do.call(rbind, norm_Vi)) # Bind to big dataframe
  colnames(bfdf) <- col.n
  bfdf$var.name <- rep(row.n,dim(varImportance)[3])
  t <- bfdf %>%
    group_by(var.name) %>%
    summarise_each(funs(mean)) #add means 
  
  return(as.data.frame(t))
  
}

normVar <- function(x,...){
  ### Internal Function for NormImp 
  out <- matrix(ncol = ncol(x),nrow = nrow(x))
  for(j in 1:ncol(x)){
    for(i in 1:nrow(x)){
      out[i,j] <- x[i,j]/sum(x[,j])
    }
  }
  return(as.data.frame(out))
}

ensVarEval <- function(ensemble.out, n){
  ###Function for retriving varriable importance scores from ensemble models
  ## Arguments:
  ## ensemble.out:  ensemble out class
  ## nrep:  number of model replicates run
  ensemble_models_names <- BIOMOD_LoadModels(ensemble.out)
  vi <- list()
  for (mod in ensemble_models_names){
    cat("\n> variables importance of ", mod)
    vi[[mod]] <- get(mod)@model_variables_importance
    #    vi <- c(vi, variables_importance(model=get(mod), data=get_formal_data(ensemble.out,'expl.var'), method="full_rand", nb_rand=n))
  }
  
  vi2 <- lapply(vi,normVar)
  vi3 <- c()
  for(i in 1:length(ensemble_models_names)){
    vi3 <- rbind(vi3,t(transmute(vi2[[i]],varI = rowMeans(vi2[[i]][,1:n]))))
  }
  rownames(vi3) <- ensemble_models_names
  colnames(vi3) <- rownames(vi[[1]])
  
  return(as.data.frame(vi3))
}

#### Batch Run Functions ####
classMod <- function(taxon.class, run.id, nrep, dbl = T, band = 4){
  ## Function for running the bioMod function with snowfall cluster sppedup
  # Arguments
  # taxon.class <- 3 letter sub set that corisponds with the breeding class
    # options are "ptr", "mol", "mic"
  # run.id <- character string used for naming the porject run
  # nrep <- number of repitions to be preformed. 
  # dbl <- Logical for if the double month varriants should be used. 
  
  #### Double month functions ####
  if(dbl){
    source("R/dblMonthFuns.R")}
  
  #### Biomod modeling function ####
  source("R/bioMod_fun.R")
  
  #### Set Up Taxon Class ####
  if(dbl){
    txc.occ <- occLoad2(taxon.class)  
  } else{
    occLoad(taxon.class)
    txc.occ <- get(paste0(taxon.class,".occ"), envir = .GlobalEnv)
  }
  
  txc.l <- listGen(txc.occ)
  
  #### Cluster Setup ####
  library(snowfall)
  sfInit(parallel=TRUE, cpus = 3)
  sfLibrary('biomod2', character.only=TRUE)
  sfLibrary('dplyr', character.only=TRUE)
  sfLibrary('RStoolbox', character.only=TRUE)
  sfLibrary('raster', character.only = TRUE)
  sfExportAll()
  
  #### bioMod Function call ####
  sfLapply(txc.l, bioMod, 
           occ.db = txc.occ,
           nrep = nrep,
           run.id = run.id,
           dbl = dbl)
  
  #### End Cluster ####
  sfStop(nostop = FALSE)
  
}
#### Post Processing ####
resRasterLoad <- function(taxon.class, run.id, dbl, path.to.dir, band = 4){
  ##Function of loading the raster results from batch runs of bioMod
  ##Arguments:
  # taxon.class <- 3 letter sub set that corisponds with the breeding class
    # options are "ptr", "mol", "mic"
  # run.id <- character string used for naming the porject run
  # dbl <- logical. If the double month modeling took place
  # path.to.dir <- path argument to the destination 
  # band <- which band to read in in refernece to the .gri file returned.
  # options are 1: EMmeanByROC, 2: EMmedianByROC, 3: EMcaByROC, 4: EMwmeanByROC
  require(raster)
  if(dbl){
    source("R/dblMonthFuns.R")
    txc.l <- listGen(occLoad2(taxon.class))
  } else { 
    occLoad(taxon.class)
    txc.l <- listGen(get(paste0(taxon.class,".occ"), envir = .GlobalEnv))}
  f.path <- file.path(path.to.dir,run.id,txc.l,txc.l,
                      paste0("proj_",txc.l,"ENSEMBLE_PROJ"),
                      paste0("proj_",txc.l,"ENSEMBLE_PROJ_",txc.l,
                             "_ensemble.gri"))
  tf.path <- f.path[which(file.exists(f.path))] #scrape out ones that DNE
  r.stk <- do.call(stack,lapply(X = tf.path,FUN = raster, band = band))
  names(r.stk) <- txc.l[which(file.exists(f.path))]
  
  return(r.stk)
}

flickerPlot <- function(res.stk, dbl, births = F, virus = F,
                        source.path = data.source, path.out = NULL){
  ##Function to plot bioMod results for .gif assembly 
  ##Arguments:
  # res.stk <- result from resRasterLoad
  # taxon.class <- 3 letter sub set that corisponds with the breeding class
    # options are "ptr", "mol", "mic"
  # dbl <- logical. If the double month modeling took place
  # births <- Logical if the birth data should be included on the map surface
  # virus <- Logical if virus data should be included on the map surface
  # Source.path <- path argument to the location of the dataframe, polygons, and
  # cropping masks needed for the figures
  # path.out <- path argument for if and where the .gif should be saved
  require(ggplot2);require(rgdal)
  
  #### Set Up ####
  ## accessory layers
  afr.poly <- readOGR(dsn = file.path(source.path, "Africa"),
                      layer = "AfricanCountires")
  rf.poly <- rasterToPolygons(raster(file.path(source.path, "cropMask.tif")),
                              fun = function(x){x==1}, dissolve = T)
  
  ## dataframe for plotting
  res.df <- data.frame(rasterToPoints(res.stk))
  
  ## name vector to itterate across 
  occ.l <- names(res.stk)
  
  #### function for ggplots ####
  gfun <- function(occ.l){
    
    #### Main Plot ####
    ## DataFrame
    t <- cbind(res.df[,1:2], res.df[,occ.l])
    colnames(t) <- c("long","lat","score")
    
    g.plot <- ggplot(t) +
      
      #create african continent background
      geom_polygon(data = fortify(afr.poly),
                   aes(long, lat, group = group), 
                   colour = "grey20",
                   alpha = .25) +
      aes(x=long, y=lat) +
      scale_fill_gradient(low = "yellow", high = "red4",
                          limits = c(0,1000))+
      geom_raster(aes(fill = score), interpolate = T)+
      
      #add area modeled
      geom_polygon(data = fortify(rf.poly),
                   aes(long, lat, group = group),
                   colour = "black", 
                   fill = NA) +
      
      #create african continent background
      geom_polygon(data = fortify(afr.poly),
                   aes(long, lat, group = group), 
                   colour = "grey20",
                   fill = NA,
                   alpha = .2) +
      ggtitle(paste0("Modeled Births ", occ.l)) +
      coord_fixed()
    
    #### Theme settings ####
    bkg <- theme(
      panel.background = element_rect(fill = "lightblue",
                                      colour = "lightblue",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = "white"),
      plot.title = element_text(hjust = 0.5))
    
    #### If options ####
      #### No births or virus ####
    if(births == F & virus == F){
      out <-  g.plot + bkg
    }
    
      #### births ####
    if(births){
      if(dbl){
        breeding.db <- occLoad2(substring(occ.l[[1]], 1, 3))
      } else {## Get brith database (should be loaded in resRasterLoad)
        breeding.db <- get(paste0(substring(occ.l[[1]], 1, 3),".occ"),
          envir = .GlobalEnv)}
      q <- subset(breeding.db[,1:2], breeding.db[,occ.l] == 1)
      q$Births <- as.factor(1)
      
      bat.births <- geom_point(data = fortify(q), aes(x = x, y = y, group = Births, shape = Births),
                               colour = "green2", 
                               size = 2, 
                               alpha = .2)
      out <- g.plot + bat.births + bkg
    }
    
      #### virus ####
    if(virus){
      vir <- read.csv(file.path(clean.dir, "humOB.PPM.csv"))
      if(dbl){
        
      } else{
        ob <- paste0("OB", substring(occ.l,4,4))
        q <- subset(vir[,2:3], vir[,ob] > 0)
        q$HumOB <- as.factor(1)
        
      }
    }
    return(out)
  }
  
  ##Apply gfun() across the occ.l
  out.l <- lapply(occ.l, gfun)
  
  #### write conditions ####
  if(!is.null(path.out)){
    ## Create path.out if it doesn't exist
    if(!dir.exists(path.out)){
      dir.create(path = path.out, recursive = T, showWarnings = F)}
    ## save as .png for ,git assembly
    #Save as .png for .gif assembly
    for(i in 1:length(out.l)){
      png(file = file.path(path.out,paste0(occ.l[[i]],".",births,virus,".png")),
          bg = "transparent", width = 680, height = 750)
      print(out.l[[i]])
      dev.off()
    }
  }
  
  return(out.l)
}

gifMaster <- function(taxon.class, path.to.sng = NULL, path.to.dbl = NULL, fps, master = T,
                      path.out = NULL){
  ## Function for creating .gif giles of flickerplots including those for 
  ## all model results together (both double and single months)
  ##Arguments:
  # taxon.class <- 3 letter sub set that corisponds with the breeding class
  # options are "ptr", "mol", "mic"
  # fps <- frames per second for gif rate. 
  # master <- Logical. if TRUE both paths need to be filled
  # path.out <- path to save the image at
  require(magick); require(gtools)
  
  
  if(master){
    v.l <- c()
    for(i in 1:12){
      pat1 <- paste0(taxon.class, i, "[T OR \\.]")
      if(length(list.files(path.to.sng, pattern = pat1)) ==1 ){
        v.l <- append(v.l, image_read(file.path(path.to.sng,
                                      list.files(path.to.sng, pattern = pat1))))}
      pat2 <- paste0(taxon.class, i,".", taxon.class, i+1)
      if(length(list.files(path.to.dbl, pattern = pat2)) ==1 ){
        v.l <- append(v.l, image_read(file.path(path.to.dbl,
                                      list.files(path.to.dbl, pattern = pat2))))}
    }
    gif <- image_animate(v.l, fps = fps)
  } else{
    ifelse(!is.null(path.to.sng),
           path.item <- path.to.sng,
           path.item <- path.to.dbl)
    p.l <- mixedsort(list.files(path.item))
    i.v <- c()
    for(i in 1:length(p.l)){
      i.v <- append(i.v, image_read(file.path(path.item, p.l[[i]])))
    }
    gif <- image_animate(i.v, fps = fps)
  }
  
  if(!is.null(path.out)){
    ## Create path.out if it doesn't exist
    if(!dir.exists(path.out)){
      dir.create(path = path.out, recursive = T, showWarnings = F)}
    if(master){
      bogie <- "MASTER"
    } else{
      if(!is.null(path.to.sng) & is.null(path.to.dbl)){
        bogie <- "_SNG"
      }
      if(!is.null(path.to.dbl) & is.null(path.to.sng)){
        bogie <- "_DBL"
      }
    }
    image_write(gif, path = file.path(path.out, paste0(taxon.class,bogie, ".gif")),
                format = "gif", quality = 100)
  }
  return(gif)
}

#### DataFrame manipulation tools ####
wrap <- function (x, n = 1L, order_by = NULL, ...){
  if (!is.null(order_by)) {
    return(with_order(order_by, wrap, x, n = n))
  }
  if (length(n) != 1 || !is.numeric(n) || n < 0) {
    dplyr:::bad_args("n", "must be a nonnegative integer scalar, ", 
                     "not {rlang::type_of(n)} of length {length(n)}")
  }
  if (n == 0) 
    return(x)
  xlen <- length(x)
  n <- n %% xlen
  out <- x[c(seq_len(n)+xlen-n, seq(xlen-n))]
  attributes(out) <- attributes(x)
  out
}
