#### Force of birthing rasters ####
source("R/helperFunctions.R"); source("R/dblMonthFuns.R")
# #### Breeding probability rasters without imputation ####
ptr.sng.raw <- resRasterLoad("ptr", "SNG_3",F,  mod.out.dir)
mic.sng.raw <- resRasterLoad("mic", "SNG_3",F,  mod.out.dir)
mol.sng.raw <- resRasterLoad("mol", "SNG_3",F,  mod.out.dir)

ptr.dbl.raw <- resRasterLoad("ptr", "DBL_3",T,  mod.out.dir)
mic.dbl.raw <- resRasterLoad("mic", "DBL_3",T,  mod.out.dir)
mol.dbl.raw <- resRasterLoad("mol", "DBL_3",T,  mod.out.dir)

#### Adding 0 where necessary ####
ptr.sng.imp <- resRasterLoad("ptr", "SNG_3",F,  mod.out.dir)
mic.sng.imp <- resRasterLoad("mic", "SNG_3",F,  mod.out.dir)
mol.sng.imp <- resRasterLoad("mol", "SNG_3",F,  mod.out.dir)

ptr.dbl.imp <- resRasterLoad("ptr", "DBL_3",T,  mod.out.dir)
mic.dbl.imp <- resRasterLoad("mic", "DBL_3",T,  mod.out.dir)
mol.dbl.imp <- resRasterLoad("mol", "DBL_3",T,  mod.out.dir)
#
## iterates through all classes, checks which ones exist, and deals with them appropriatly where they do not
tax <- c("ptr", "mol", "mic")
rf <- raster(file.path(data.source, "cropMask.tif"))
q <- reclassify(rf, rcl = c(0,1,0))
lapply(tax, occLoad)
## fill 0 ##
for(i in 1:length(tax)){
  #Single
  df <- get(paste0(tax[[i]],".occ"))
  l <- names(df[,3:ncol(df)]) #names of what should be there
  z <- which(l %!in% names(get(paste0(tax[[i]], ".sng.imp")))) #which aren't modeled
  if(length(z) > 0){
    z.empt <- do.call(stack, replicate(length(z),q))
    names(z.empt) <- l[z]
    assign(paste0(tax[[i]], ".sng.imp"), stack(get(paste0(tax[[i]], ".sng.imp")), z.empt))
  }
  #Double
  df <- occLoad2(paste0(tax[[i]]))
  l <- c()#names of what should be there
  for(g in 1:12){
    v <- 1:12
    j <- ((v + g) - 2) %% length(v)+1
    k <- j[1:2]
    m1 <- paste0(tax[[i]],k[[1]])
    m2 <- paste0(tax[[i]],k[[2]])
    l[[g]] <- paste0(m1,".",m2)
  }
  z <- which(l %!in% names(get(paste0(tax[[i]], ".dbl.imp"))))
  if(length(z) > 0){
    z.empt <- do.call(stack, replicate(length(z),q))
    names(z.empt) <- l[z]
    assign(paste0(tax[[i]], ".dbl.imp"), stack(get(paste0(tax[[i]], ".dbl.imp")), z.empt))
  }
}

div.stk <- do.call(stack,
                   lapply(file.path(clean.dir, list.files(clean.dir, pattern = "*.div.tif")[c(2,3,5)]),
                           raster))


 ## Changes names of stack layers
   ##Note: Sadly does not work with lapply of with for loops
 nmr <- function(x){
   ## Function for renaming the raster stacks for further processing
   lz <- list()
     for(i in 1:nlayers(x)){
     lz[[i]] <- paste(deparse(substitute(x)),substr(names(x[[i]]), 4,nchar(names(x[[i]]))), collapse = "_" )
   }
   names(x) <- unlist(lz)
   return(x)
 }

 ptr.sng.raw <- nmr(ptr.sng.raw)
 ptr.sng.imp <- nmr(ptr.sng.imp)
 ptr.dbl.raw <- nmr(ptr.dbl.raw)
 ptr.dbl.imp <- nmr(ptr.dbl.imp)

 mic.sng.raw <- nmr(mic.sng.raw)
 mic.sng.imp <- nmr(mic.sng.imp)
 mic.dbl.raw <- nmr(mic.dbl.raw)
 mic.dbl.imp <- nmr(mic.dbl.imp)

 mol.sng.raw <- nmr(mol.sng.raw)
 mol.sng.imp <- nmr(mol.sng.imp)
 mol.dbl.raw <- nmr(mol.dbl.raw)
 mol.dbl.imp <- nmr(mol.dbl.imp)

 ### list stacks ###
 p.stks <- c(ptr.sng.raw, ptr.sng.imp, ptr.dbl.raw, ptr.dbl.imp)
 mi.stks <- c(mic.sng.raw, mic.sng.imp, mic.dbl.raw, mic.dbl.imp)
 mo.stks <- c(mol.sng.raw, mol.sng.imp, mol.dbl.raw, mol.dbl.imp)

 br.force <- function(x, y, write = FALSE){
   ## function for adding the force of birthing to all!
   ## N.B. This is not the same log +1 force used in the model but raw
   lz.out <- list()
   for( i in 1:nlayers(x)){
     lz.out[[i]] <- x[[i]]/1000 * y
   }
   stk.out <- stack(lz.out)
   names(stk.out) <- names(x)

   if(write == T){
     writeRaster(stk.out,
                 filename = file.path(clean.dir, "BirthForce", "BR"),
                 format = "GTiff",
                 bylayer = T,
                 suffix = "names",
                 overwrite = T)
   }

   return(stk.out)
 }

 ## add force of birthing
 p.stks.br <- lapply(p.stks, br.force, y = div.stk$ptr.div, write = T)
 mi.stks.br <- lapply(mi.stks, br.force, y = div.stk$mic.div, write = T)
 mo.stks.br <- lapply(mo.stks, br.force, y = div.stk$mol.div, write = T)

## check to make sure the length of above is correct ##
# sum(unlist(lapply(c(p.stks,mo.stks, mi.stks),function(x)(sum(nlayers(x))))))
# length(list.files(file.path(clean.dir,"BirthForce")))
#### Checkpoint 1, stacks above were writen out and saved ####
