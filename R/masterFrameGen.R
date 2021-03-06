#### Creating the longMaster dataframe, take 2 ####
## run to create the dataframe required for regression analysis

source("R/helperFunctions.R"); source("R/dblMonthFuns.R")

## Need to create long format dataframe for regression analysis of ebola outbreak events 
source("R/helperFunctions.R"); source("R/dblMonthFuns.R")
#### include 0 where there is not any data
#### OutBreak Events ####
# outbreak events have previously been split into human/and animal events, and seperateded across months
library(raster);library(gtools)
# human outbreak stack
hum.stk0 <- do.call(stack, lapply(
  file.path(clean.dir.nov,mixedsort(list.files(file.path(clean.dir.nov), pattern = "OB_hum0*")))
  ,raster))

# +2 human outbreak stack
hum.stk0_2 <- do.call(stack, lapply(
  file.path(clean.dir.nov,mixedsort(list.files(file.path(clean.dir.nov), pattern = "OB_2_hum0_*")))
  ,raster))

# animal outbreak stack
ann.stk0 <- do.call(stack, lapply(
  file.path(clean.dir.nov,mixedsort(list.files(file.path(clean.dir.nov), pattern = "OB_ann0_*")))
  ,raster))

# host detection locations (bats with detected virus)
hdl.stk0 <- do.call(stack, lapply(
  file.path(clean.dir.nov,mixedsort(list.files(file.path(clean.dir.nov), pattern = "OB_hdl0_*")))
  ,raster))

#### EMN breeding probability rasters with 0 back imputation #####
ptr.dbl.imp <- calc(resRasterLoad("ptr", "DBL_3",T,  mod.out.dir),  function(x){x/1000})
mic.dbl.imp <- calc(resRasterLoad("mic", "DBL_3",T,  mod.out.dir), function(x){x/1000})
mol.dbl.imp <- calc(resRasterLoad("mol", "DBL_3",T,  mod.out.dir), function(x){x/1000})


### Taxon stuff
tax <- c("ptr", "mol", "mic")
rf <- raster(file.path(data.source, "cropMask.tif"))
q <- reclassify(rf, rcl = c(0,1,0))
lapply(tax, occLoad)

## fill 0  step##
for(i in 1:length(tax)){
  ## fill in the backrougn mask where model resuts were not around
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
                   lapply(file.path(clean.dir, list.files(clean.dir, pattern = "*.div.tif")[2:5]),
                          raster))
frag <- raster(file.path(clean.dir, "fragIndex.tif"))
sttc.stk <- do.call(stack, c(pop.den, land.cover, div.stk, frag))

#### Creating the two month windows ####
long.list <- list()
blank <- rf
values(blank) <- NA
for(i in 1:12){
  #Outbreak Stack
  ob.stk <- stack(hum.stk0[[paste0("OB_hum0_",i)]],
                  hum.stk0_2[[paste0("OB_2_hum0_",i)]],
                  ann.stk0[[paste0("OB_ann0_",i)]],
                  hdl.stk0[[paste0("OB_hdl0_",i)]]) 
  names(ob.stk) <- c("OB_hum_imp","OB_hum_imp_2", "OB_ann_imp", "hdl")
  #dbl month of breeding impute
  br.dbl.imp <- list()
  for(j in 1:3){
    tax <- c(ptr.dbl.imp, mic.dbl.imp, mol.dbl.imp)
    tax.name <- c("ptr", "mic", "mol")
    # Need to use the current month, plus the previous month here.
    current_month <- i
    previous_month <- if (i == 1) 12 else i-1
    dbl_name <- paste0(tax.name[[j]],previous_month,".",tax.name[[j]],current_month)
    if(dbl_name %in% names(tax[[j]])){
      br.dbl.imp[[j]] <- tax[[j]][[dbl_name]]
    } else {
      br.dbl.imp[[j]] <- blank
    }
  }
  dbl.stk.imp <- do.call(stack, br.dbl.imp)
  names(dbl.stk.imp) <- c("ptr_dbl_imp", "mic_dbl_imp", "mol_dbl_imp")
  #Bring it all together
  month.stk <- do.call(stack, c(ob.stk, dbl.stk.imp, sttc.stk))
  month.df <- as.data.frame(month.stk)
  month.df$cell <- paste0("c",seq(1:nrow(month.df)))
  month.df$month <- i
  long.list[[i]] <- month.df
}

#### COmbine and create loged terms
library(data.table);library(tidyverse)
long.table <- as.data.table(do.call(rbind, long.list))
xy <- as.data.table(xyFromCell(blank, seq(1:ncell(blank))))
xy$cell <- paste0("c", seq(1:ncell(blank)))
long.table <- left_join(long.table, xy, "cell") %>%
  ## create additional variables
  mutate(mlptr_prob = log(ifelse(y > -10 & month >7, # modified Ptr layer
                               (5/9*ptr_dbl_imp + 1),
                               ptr_dbl_imp + 1)),
         mptr_prob = ifelse(y > -18 & month >7, # modified Ptr layer
                               (5/9*ptr_dbl_imp + 1),
                               ptr_dbl_imp + 1),
         BB.cond = (1-((1-ptr_dbl_imp)*(1-mic_dbl_imp)*(1-mol_dbl_imp))), #Conditions
         lBB.cond = log(BB.cond + 1),
         mBB.cond = (1-((1-mptr_prob)*(1-mic_dbl_imp)*(1-mol_dbl_imp))),
         mlBB.cond = log(mBB.cond + 1),
         BatDiv = (ptr.div + mic.div + mol.div), #Totaled bat Diversity
         lBatDiv = log(BatDiv +1),
         BB.condDIV = BB.cond * BatDiv,
         mBB.condDIV = mBB.cond * BatDiv,
         lBB.condDIV = log(BB.condDIV + 1),
         mlBB.condDIV = log(mBB.condDIV + 1),
         lptr_BR = log(ptr_dbl_imp * ptr.div + 1), #force of birthing
         mlptr_BR = log(mptr_prob * ptr.div + 1),
         lmic_BR = log(mic_dbl_imp * mic.div + 1),
         lmol_BR = log(mol_dbl_imp * mol.div + 1),
         lptr_prob = log(ptr_dbl_imp + 1),
         lmic_prob = log(mic_dbl_imp + 1),
         lmol_prob = log(mol_dbl_imp + 1),
         lPop = log(popDen + 1),
         lFrag = log(fragIndex + 1),
         lnBm.div = log(nBm.div +1))
# Warning messages:
#   1: In log(BB.cond + 1) : NaNs produced
#   2: In log(BB.condDIV + 1) : NaNs produced

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

## create list of items that need to be lagged and then use a function to apply accross
## log taxon names

wrapperTrapper <- function(x, df = long.table, depth){
  ## function to allow automated wrapping to multipul columns at the same time 
  ## to a consistant interval
  x.enquo <- enquo(x)
  col.out.names <- list()
  lag.columns <- list()
  q <- 1
  for(i in 1:depth){
    new.name <- paste(quo_name(x.enquo), i, sep="_")
    bz <- df %>%
      dplyr::group_by(cell) %>%
      
      dplyr::mutate_(.dots = setNames(list(lazyeval::interp(~wrap(x = nw,
                                                                  n = p,
                                                                  order_by = month),
                                                            nw=as.name(x),
                                                            p=i,
                                                            month = as.name("month"))),
                                      new.name)) %>%
      ungroup %>%
      dplyr::select(!!new.name)
    lag.columns[[q]] <- bz
    col.out.names[[q]] <- new.name
    q <- q+1
    cat(new.name,"\n")
  }
  out <- do.call(bind_cols, lag.columns)
  names(out) <- col.out.names
  return(out)
  
}

## THings that need to be lagged
## birth pulses all styles
## animal outbreaks

ltax <- paste0("l",tax.name)
birthPulses <- names(long.table)[unlist(lapply(ltax, grep, x = names(long.table)))]
conditionals <- names(long.table)[grep("lBB", x = names(long.table))]


## lag to six months 
lag6 <- c(birthPulses, conditionals)

lag6.df <- lapply(lag6,
                  wrapperTrapper,
                  depth = 6)

## lag to one month 
lag1 <- wrapperTrapper("OB_ann_imp",
                       df = long.table,
                       depth = 1)
longMasterFull <- do.call(bind_cols, c(long.table, lag6.df, lag1))
longMasterFull <- longMasterFull %>%
  mutate(lptrDiv = log(ptr.div + 1), ## add logged split Div so it doesnt lag
         lmicDiv = log(mic.div + 1),
         lmolDiv = log(mol.div + 1))
fwrite(longMasterFull, 
       file = file.path(clean.dir.nov,"fullLongMaster_Complete.csv"))
