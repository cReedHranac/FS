## Creating the Super frame 
clean.env <- list(ls())

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

#### Breeding probability rasters without imputation #### 
ptr.sng.raw <- resRasterLoad("ptr", "SNG_2",F,  mod.out.dir)
mic.sng.raw <- resRasterLoad("mic", "SNG_2",F,  mod.out.dir)
mol.sng.raw <- resRasterLoad("mol", "SNG_2",F,  mod.out.dir)

ptr.dbl.raw <- resRasterLoad("ptr", "DBL_2",T,  mod.out.dir)
mic.dbl.raw <- resRasterLoad("mic", "DBL_2",T,  mod.out.dir)
mol.dbl.raw <- resRasterLoad("mol", "DBL_2",T,  mod.out.dir)

#### Adding 0 where necessary ####
ptr.sng.imp <- resRasterLoad("ptr", "SNG_2",F,  mod.out.dir)
mic.sng.imp <- resRasterLoad("mic", "SNG_2",F,  mod.out.dir)
mol.sng.imp <- resRasterLoad("mol", "SNG_2",F,  mod.out.dir)

ptr.dbl.imp <- resRasterLoad("ptr", "DBL_2",T,  mod.out.dir)
mic.dbl.imp <- resRasterLoad("mic", "DBL_2",T,  mod.out.dir)
mol.dbl.imp <- resRasterLoad("mol", "DBL_2",T,  mod.out.dir)

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

#### Static Co-varriates ####
# Population density, Biodiversity, Landcover, Fragmentation
pop.den <- raster(file.path(clean.dir, "popDen.tif"))
land.cover <- raster(file.path(clean.dir, "LandCover.tif"))
div.stk <- do.call(stack, 
                   lapply(file.path(clean.dir, list.files(clean.dir, pattern = "*_sum.tif")),
                          raster))
frag <- raster(file.path(clean.dir, "fragIndex.tif"))
sttc.stk <- do.call(stack, c(pop.den, land.cover, div.stk, frag))

## longFrame ## 
long.list <- list()
blank <- rf
values(blank) <- NA
for(i in 1:12){
  #Outbreak Stack
  ob.stk <- stack(hum.stk0[[paste0("OB_hum0_",i)]], ann.stk0[[paste0("OB_ann0_",i)]],
                  hum.stk[[paste0("OB_hum",i)]], ann.stk[[paste0("OB_ann",i)]]) 
  names(ob.stk) <- c("OB_hum_imp", "OB_ann_imp", "OB_hum_raw", "OB_ann_raw")
  #Breeding Stacks
  ####Raw####
    #sng
  br.sng.raw <- list()
  for(j in 1:3){
    tax <- c(ptr.sng.raw, mic.sng.raw, mol.sng.raw)
    tax.name <- c("ptr", "mic", "mol")
    if(paste0(tax.name[[j]],i) %in% names(tax[[j]])){
      br.sng.raw[[j]] <- tax[[j]][[paste0(tax.name[[j]],i)]]
    } else {
      br.sng.raw[[j]] <- blank
    }
  }
  sng.raw.stk <- do.call(stack, br.sng.raw)
  names(sng.raw.stk) <- c("ptr_sng_raw", "mic_sng_raw", "mol_sng_raw")
    #dbl
  br.dbl.raw <- list()
  for(j in 1:3){
    tax <- c(ptr.dbl.raw, mic.dbl.raw, mol.dbl.raw)
    tax.name <- c("ptr", "mic", "mol")
    v <- 1:12
    k <- ((v + i) - 2) %% length(v)+1
    l <- k[2]
    if(paste0(tax.name[[j]],i,".",tax.name[[j]],l) %in% names(tax[[j]])){
      br.dbl.raw[[j]] <- tax[[j]][[paste0(tax.name[[j]],i,".",tax.name[[j]],l)]]
    } else {
      br.dbl.raw[[j]] <- blank
    }
  }
  dbl.stk.raw <- do.call(stack, br.dbl.raw)
  names(dbl.stk.raw) <- c("ptr_dbl_raw", "mic_dbl_raw", "mol_dbl_raw")
  ####INPT####
  #sng
  br.sng.imp <- list()
  for(j in 1:3){
    tax <- c(ptr.sng.imp, mic.sng.imp, mol.sng.imp)
    tax.name <- c("ptr", "mic", "mol")
    if(paste0(tax.name[[j]],i) %in% names(tax[[j]])){
      br.sng.imp[[j]] <- tax[[j]][[paste0(tax.name[[j]],i)]]
    } else {
      br.sng.imp[[j]] <- blank
    }
  }
  sng.imp.stk <- do.call(stack, br.sng.imp)
  names(sng.imp.stk) <- c("ptr_sng_imp", "mic_sng_imp", "mol_sng_imp")
  #dbl
  br.dbl.imp <- list()
  for(j in 1:3){
    tax <- c(ptr.dbl.imp, mic.dbl.imp, mol.dbl.imp)
    tax.name <- c("ptr", "mic", "mol")
    v <- 1:12
    k <- ((v + i) - 2) %% length(v)+1
    l <- k[2]
    if(paste0(tax.name[[j]],i,".",tax.name[[j]],l) %in% names(tax[[j]])){
      br.dbl.imp[[j]] <- tax[[j]][[paste0(tax.name[[j]],i,".",tax.name[[j]],l)]]
    } else {
      br.dbl.imp[[j]] <- blank
    }
  }
  dbl.stk.imp <- do.call(stack, br.dbl.imp)
  names(dbl.stk.imp) <- c("ptr_dbl_imp", "mic_dbl_imp", "mol_dbl_imp")
  #Bring it all together
  month.stk <- do.call(stack, c(ob.stk, sng.raw.stk, dbl.stk.raw, sng.imp.stk, dbl.stk.imp, sttc.stk))
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



#### Add force of breeding ####
tax <- c("ptr", "mic", "mol") #taxonomic group
tx <- c("Mega", "Micro", "Molo")
grp <- c("sng", "dbl") #temporal grouping
hndl <- c("raw", "imp") #handle as raw or imputed

item <- list()
q <- 1
for(i in 1:length(tax)){ #tax
  for(j in 1:length(grp)){ #group
    for(k in 1:length(hndl)){ #handle
      nm <- paste(tax[[i]], grp[[j]], hndl[[k]], "BR",sep="_")
      breed.prob <- paste(tax[[i]], grp[[j]], hndl[[k]], sep="_")
      class.div <- paste0(tx[[i]],"_sum")
      item[[q]] <- long.table %>% 
        mutate_(.dots = setNames(list(lazyeval::interp(~log(breed.prob * class.div + 1),
                                                       breed.prob = as.name(breed.prob),
                                                       class.div = as.name(class.div))), nm)) %>%
        dplyr::select(!!nm)
      q <- q+1
    }
  } 
}
f.br <- do.call(cbind, item)
long.table.br <- as.data.table(cbind(long.table, f.br))
# Warning message:
#   In data.row.names(row.names, rowsi, i) :
#   some row.names duplicated: 3,4,5,
#### Add lags ####
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

q <- 1
item <- list()
## add lag columns ##
for(i in 1:length(tax)){ #tax
  for(j in 1:length(grp)){ #group
    for(k in 1:length(hndl)){ #handle
      for(l in 1:6){ # n lag
        nm <- paste(tax[[i]], grp[[j]], hndl[[k]], "BR", l ,sep="_")
        nw <- paste(tax[[i]], grp[[j]], hndl[[k]], "BR",sep="_")
        bz <- long.table.br %>% dplyr::group_by(cell) %>%
          dplyr::mutate_(.dots = setNames(list(lazyeval::interp(~wrap(x = nw, n = p, order_by = month),
                                                                nw=as.name(nw),
                                                                p=l,
                                                                month = as.name("month"))),nm)) %>%
          ungroup %>%
          dplyr::select(!!nm)
        item[[q]] <- bz
        q <- q+1
        cat(nm,"\n")
      }
    }
  } 
}

lagz <- as.data.table(do.call(cbind, item))
long.table.full <- long.table.br %>%
  bind_cols(lagz) %>%
  mutate(logPop = log(popDen + 1),
         NB_lDiv = log((Mam_sum - (Mega_sum + Molo_sum + Micro_sum))+1))

fwrite(long.table.full, file.path(clean.dir, "longMaster.csv"))
#Clean env
rm(list = setdiff(clean.env,ls()))
