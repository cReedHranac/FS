## NormSmash
  ## Script for normalizing temporal varriables

#Rescale function
rescale.raster <- function(data0, data1=NULL){
  if(missing(data1)){
    raster <- (data0-cellStats(data0,'mean'))/cellStats(data0,'sd')
  }else{
    raster <- (data1-cellStats(data0,'mean'))/cellStats(data0,'sd')
  }
  return(raster)
}

#Writer function
normSmash <- function(x){
  require(raster); require(rgdal)
  ## A function to read in, normalize and write out temporal co-varriates
  
  #Read in 
  in.stk <- stack(lapply(file.path(clean.dir,
                                   list.files(file.path(clean.dir),pattern = c(x, "*"))),
                         raster))
  
  #Normalize
  norm.stk <- rescale.raster(in.stk)
  
  #Write out 
  writeRaster(norm.stk, filename = paste0(file.path(norm.dir),"/norm.tif"),
              format = "GTiff", bylayer = T, suffix = "names")
}

## Create vector of temporal co-variates
tempo <- c("evi", "pet", "prec", "tmean")

## Check if they already exist and create with normSmash if not

lapply(tempo, function(x){
  if(length(file.exists(file.path(norm.dir,
    list.files(file.path(norm.dir), pattern = c(x,"*"))))) == 12){
    cat("Files for", x, "already exist\n")
    }else {cat("Files for", x, "need to be generated\n")
      normSmash(x)}})
         

