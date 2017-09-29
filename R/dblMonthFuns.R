#################################################
##### DoubleMonth Functions for FiloSpatial #####
##### Re-Script Sept 2017                   #####
#################################################

#### Pre-Processing Functions ####
occLoad2 <- function(taxon.class){
  ### Function to combine months into 2 month clusters of with sufficient data
  ### for modeling 
  ## Agrugments:
  ## taxon.class <- 3 letter sub set that corisponds with the breeding class
  # options built in are "ptr", "mol", "mic" (and "all" not currently supported)
  occLoad(taxon.class = taxon.class)
  occ.df <- get(paste0(taxon.class,".occ"), envir = .GlobalEnv)
  
  out.xy <- occ.df[,1:2]
  out.frame <- out.xy
  for(i in 1:12){
    v <- 1:12
    j <- ((v + i) - 2) %% length(v)+1
    k <- j[1:2]
    m1 <- paste0(taxon.class,k[[1]])
    m2 <- paste0(taxon.class,k[[2]])
    if(sum(occ.df[,m1],occ.df[,m2]) >= 5){
      months.bi <- rep(NA, nrow(occ.df))
      for(d in 1:nrow(occ.df)){
        ifelse(occ.df[d,m1]+occ.df[d,m2] >=1 , months.bi[d] <- 1, months.bi[d] <- 0)
      }
      out.frame <- cbind(out.frame, months.bi)
      names(out.frame) <- c(head(names(out.frame),-1), paste0(m1,".",m2))
    }
  }
  return(out.frame)
}

stkRevolver2 <- function(month, path.to.dir = norm.dir){
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
    j <- i[c(1:4,11:12)] ## Amendment for 2 data zone
    
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

