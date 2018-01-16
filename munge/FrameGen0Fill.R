## Creating the Super frame 
## Need to create long format dataframe for regression analysis of ebola outbreak events 
source("R/helperFunctions.R"); source("R/dblMonthFuns.R")
#### include 0 where there is not any data

#### OutBreak Events ####
# outbreak events have previously been split into human/and animal events, and seperateded across months
library(raster);library(gtools)
# human outbreak stack
hum.stk0 <- do.call(stack, lapply(
  file.path(clean.dir,mixedsort(list.files(file.path(clean.dir), pattern = "OB_hum0*")))
  ,raster))
hum.stk <- do.call(stack, lapply(
  file.path(clean.dir,mixedsort(list.files(file.path(clean.dir), pattern = "OB_hum*")))
  ,raster))

# animal outbreak stack
ann.stk0 <- do.call(stack, lapply(
  file.path(clean.dir,mixedsort(list.files(file.path(clean.dir), pattern = "OB_ann0*")))
  ,raster))
ann.stk <- do.call(stack, lapply(
  file.path(clean.dir,mixedsort(list.files(file.path(clean.dir), pattern = "OB_ann*")))
  ,raster))

#### Breeding probability rasters ####
ptr.sng <- resRasterLoad("ptr", "SNG_2",F,  mod.out.dir)
mic.sng <- resRasterLoad("mic", "SNG_2",F,  mod.out.dir)
mol.sng <- resRasterLoad("mol", "SNG_2",F,  mod.out.dir)

ptr.dbl <- resRasterLoad("ptr", "DBL_2",T,  mod.out.dir)
mic.dbl <- resRasterLoad("mic", "DBL_2",T,  mod.out.dir)
mol.dbl <- resRasterLoad("mol", "DBL_2",T,  mod.out.dir)

#### Adding 0 where necessary ####
tax <- c("ptr", "mol", "mic")
rf <- raster(file.path(data.source, "cropMask.tif"))
q <- reclassify(rf, rcl = c(0,1,0))
lapply(tax, occLoad)

for(i in 1:length(tax)){
#Single
  df <- get(paste0(tax[[i]],".occ"))
  l <- names(df[,3:ncol(df)]) #names of what should be there
  z <- which(l %!in% names(get(paste0(tax[[i]], ".sng")))) #which aren't modeled
  if(length(z) > 0){
    z.empt <- do.call(stack, replicate(length(z),q))
    names(z.empt) <- l[z]
    assign(paste0(tax[[i]], ".sng"), stack(get(paste0(tax[[i]], ".sng")), z.empt))  
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
  z <- which(l %!in% names(get(paste0(tax[[i]], ".dbl"))))
  if(length(z) > 0){
    z.empt <- do.call(stack, replicate(length(z),q))
    names(z.empt) <- l[z]
    assign(paste0(tax[[i]], ".dbl"), stack(get(paste0(tax[[i]], ".dbl")), z.empt))  
  }
}

#### Static Co-varriates ####
# Population density, Biodiversity, Landcover, Fragmentation
pop.den <- raster(file.path(clean.dir, "popDen.tif"))
land.cover <- raster(file.path(clean.dir, "LandCover.tif"))
div.stk <- do.call(stack, 
                   lapply(file.path(clean.dir, list.files(clean.dir, pattern = "*_sum.tif")),
                          raster))
frag <- raster(file.path(clean.dir, "fragIndex.tif"))
sttc.stk <- do.call(stack, c(pop.den, land.cover, div.stk, frag))

#### Master Stack/DataFrame ####
# super.stk <- do.call(stack, c(hum.stk, ann.stk,
#                              ptr.sng, ptr.dbl,
#                              mic.sng, mic.dbl,
#                              mol.sng, mol.dbl,
#                              pop.den, land.cover, div.stk))
# 
# super.frame <- as.data.frame(super.stk)
# write.csv(super.frame, file.path(clean.dir,"wideTable.csv"), row.names = F)
## Not what we need

## longFrame ## 
long.list <- list()
blank <- rf
values(blank) <- NA
for(i in 1:12){
  #Outbreak Stack
  ob.stk <- stack(hum.stk0[[paste0("OB_hum0_",i)]], ann.stk0[[paste0("OB_ann0_",i)]],
                  hum.stk[[paste0("OB_hum",i)]], ann.stk[[paste0("OB_ann",i)]]) 
  names(ob.stk) <- c("OB_hum0_", "OB_ann0_", "OB_hum", "OB_ann")
  #Breeding Stacks
  #sng
  br.sng <- list()
  for(j in 1:3){
    tax <- c(ptr.sng, mic.sng, mol.sng)
    tax.name <- c("ptr", "mic", "mol")
    if(paste0(tax.name[[j]],i) %in% names(tax[[j]])){
      br.sng[[j]] <- tax[[j]][[paste0(tax.name[[j]],i)]]
    } else {
      br.sng[[j]] <- blank
    }
  }
  sng.stk <- do.call(stack, br.sng)
  names(sng.stk) <- c("ptr_sng", "mic_sng", "mol_sng")
  #dbl
  br.dbl <- list()
  for(j in 1:3){
    tax <- c(ptr.dbl, mic.dbl, mol.dbl)
    tax.name <- c("ptr", "mic", "mol")
    v <- 1:12
    k <- ((v + i) - 2) %% length(v)+1
    l <- k[2]
    if(paste0(tax.name[[j]],i,".",tax.name[[j]],l) %in% names(tax[[j]])){
      br.dbl[[j]] <- tax[[j]][[paste0(tax.name[[j]],i,".",tax.name[[j]],l)]]
    } else {
      br.dbl[[j]] <- blank
    }
  }
  dbl.stk <- do.call(stack, br.dbl)
  names(dbl.stk) <- c("ptr_dbl", "mic_dbl", "mol_dbl")
  
  #Bring it all together
  month.stk <- do.call(stack, c(ob.stk, sng.stk, dbl.stk, sttc.stk))
  month.df <- as.data.frame(month.stk)
  month.df$cell <- paste0("c",seq(1:nrow(month.df)))
  month.df$month <- i
  long.list[[i]] <- month.df
}

library(data.table);library(dplyr)
long.table <- as.data.table(do.call(rbind, long.list))
xy <- as.data.table(xyFromCell(blank, seq(1:ncell(blank))))
xy$cell <- paste0("c", seq(1:ncell(blank)))
long.table <- left_join(long.table, xy, "cell")

fwrite(long.table, file.path(clean.dir, "longTable0Fill.csv"))
