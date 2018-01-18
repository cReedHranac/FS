#### spatGLM clean with imputed data ####
source("R/helperFunctions.R")
library(data.table); library(dplyr); library(tidyr)

## Dave run this!
#clean.dir <- "."

#### Functions ####
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

pppWeights <- function(ob.col, dataFrame, rasterGrid){
  require(rlang)
  ##Funciton for creating the spatial weights to be applied for the glm
  ##Arguments:
  # ob.col <- colname to be selected from dataFrame
  # dataFrame <- the long dataFrame containing all pertinent information
  # rasterGrid <- the raster to grid the results to
  
  ####POINTS ####
  ob.quo <- enquo(ob.col)
  ob.pts <-dataFrame %>% filter(UQ(ob.quo) == 1)#select poitns
  coordinates(ob.pts) <- ~ x + y
  proj4string(ob.pts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  ob.pts <- remove.duplicates(ob.pts) #  filter unique points
  ob <- ppp(ob.pts@coords[,"x"],ob.pts@coords[,"y"],
            window = windows[[1]])
  
  #### DETERMINING WEIGHTS etc####
  library(spatstat)
  im_win <- as.im(rasterGrid)
  quad_scheme <- pixelquad(ob, im_win)
  p_for_weights <- ppm(quad_scheme)
  
  # Basically we run ppm to get the weights and stuff
  P <- union.quad(p_for_weights$Q) # this combines data and dummy points (grid)
  X <- p_for_weights$Q$data        # this is just the data
  Q <- p_for_weights$Q             # this is the quadrature object (data+dummy points etc)
  nQ <- n.quad(Q)                  # number of points
  zQ <- is.data(Q)                 # This is where we have data in P (combined data/dummy)
  wQ <- w.quad(Q)                  # This is the weights
  yQ <- numeric(nQ)                # Y - outcome variable after transformation (so it works with a GLM using quasi family stuff)
  yQ[zQ] <- 1/wQ[zQ]               # Fill in Y with 1/weights where we have data
  mQ <- marks.quad(Q)              # in case we have marks (we don't)
  zeroes <- attr(wQ, "zeroes")
  sQ <- if(is.null(zeroes)) rep.int(TRUE, nQ) else !zeroes # Checking whether we are running on a subset of the region (we're not)
  
  .mpl <- list(W      = wQ,
               Z      = zQ,
               Y      = yQ,
               MARKS  = mQ,
               SUBSET = sQ)
  glmdata <- data.frame(.mpl.W = .mpl$W, .mpl.Y = .mpl$Y, x = P$x, y = P$y)
  # Expand out so we have month as well
  all <- list()
  for (i in 1:12) {
    all[[i]] <- cbind(glmdata, month=i)
  }
  glm_all <- do.call(rbind, all)
  
  return(glm_all)
}

weightedDf <- function(ob.col, obWeights.df, dataFrame, cols){
  ##Function for created and dataframe with weighted responce
  ##Arguments:
  # ob.col <- colname to be selected from dataFrame
  # obWeights.df <- dataframe returned by pppWeights
  # dataFrame <- dataframe with covariate data (also used in pppWeights)
  # cols <- list of columns to be included from dataFrame
  gridMe <- function(x, mx2) {
    min_x <- min(x,mx2)
    dx <- diff(sort(unique(x)))
    delta_x <- min(dx[dx > 1e-10])
    round((x - min_x) / delta_x)
  }
  #### Select columns ####
  dat.sub <- dataFrame %>%
    dplyr::select(cols) #select select data
  #### handle the gridding being different ####
  obGrid <- obWeights.df %>% 
    mutate_(x_grid = ~gridMe(x, min(dataFrame[,"x"])), ##From the raster grid
            y_grid = ~gridMe(y, min(dataFrame[,"y"])))
  datGrid <- dat.sub %>% 
    mutate_(x_grid = ~gridMe(x, min(obWeights.df[,"x"])), ##From the ppp grid
            y_grid = ~gridMe(y, min(obWeights.df[,"y"]))) %>%
    dplyr::rename_(x2=~x, y2=~y)
  
  if(!any(obGrid$x_grid %in% datGrid$x_grid) ||
     !any(obGrid$y_grid %in% datGrid$y_grid)){
    stop("You messed something up with the grid stuff, better call Jonathan/n")
  }
  #### Combine dataframes ####
  frameFull <- left_join(obGrid, datGrid, by = c("x_grid", "y_grid", "month")) 
  quo.col <- enquo(ob.col)
  frame.w <- frameFull %>% 
    mutate(.mpl.Y = .mpl.Y * !!quo.col) %>% 
    replace_na(list(.mpl.Y = 0))
  return(frame.w)
}

spatGLM <- function(ob.col, coV.v, dat, rGrid = rf){
  ###Function for applying the hybrid spatGLM to data
  ###Arguments
  ### ob.col <- outbreak column from dataframe
  ### coV.v <- vector of covariate names from dat
  # Note: last 5 items should include month, ob.col, x,y, and cell
  ### dat <- dataframe containing the relivent information
  ### rGrid <- rasterGrid from pppWeights obj (defalut to rf mask)
  #### creating point object ####
  ob.col <- enquo(ob.col)
  ppp.ob <- pppWeights(ob.col = UQ(ob.col),
                       dataFrame = dat,
                       rasterGrid = rGrid)
  #### Weighted dataframe/ weights vector ####
  W.df <- weightedDf(ob.col = UQ(ob.col),
                     obWeights.df = ppp.ob,
                     dataFrame = dat,
                     cols = coV.v)
  W.v <- W.df$.mpl.W #weights vector
  #### Model ####
  form <- as.formula(paste(".mpl.Y ","~",
                           paste(coV.v[1:(length(coV.v)-5)],collapse = "+")))
  
  mod <- glm(form,
             family=quasi(link="log", variance="mu"),
             weights=W.v,
             data=W.df)
  #### Predict ####
  pred <- predict.glm(mod,newdata = W.df, na.action = na.pass)
  pred.df <- cbind(pred, W.df)
  
  rast <- list()
  for( i in 1:12){
    empty.raster <- rGrid
    empty.raster[] <- NA_real_
    cell.dance <- pred.df %>%
      filter(month == i)
    empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$pred
    names(empty.raster) <- paste0("pred_",i)
    rast[[i]] <- empty.raster
  }
  
  items.out <- list(mod, pred.df, rast)
}
#### Data ####
dat <- tbl_df(fread(file.path(clean.dir, "longTable0Fill.csv")))
# library(skimr)
#  skim(dat) 
  ## Add force of breeding ##
 tax <- c("ptr", "mic", "mol") #taxonomic group
 tx <- c("Mega", "Micro", "Molo")
 grp <- c("sng", "dbl") #temporal grouping
 hndl <- c("raw", "imp") #handle as raw or imputed
 
item <- list()
q <- 1
 for(i in 1:length(tax)){ #tax
   for(j in 1:length(grp)){ #group
     for(k in 1:length(hndl)){ #handle
       nm <- paste(tax[[i]], grp[[j]], hndl[[k]], "BR",sep=".")
       z <- ifelse(k == 1, log(dat[,paste(tax[[i]], grp[[j]],sep="_")] * dat[,paste0(tx[[i]],"_sum")] +1),
                               log(dat[,paste0(tax[[i]],"_", grp[[j]],"I")] * dat[,paste0(tx[[i]],"_sum")] +1))
       item[[q]] <- z ; names(item[[q]])  <- nm
       q <- q+1
     }
   } 
 }
 
 mutator <- function(tax, df){
   tax.q <- enquo(tax)
   
 }
 
 dat.br <- dat %>%
  mutate(ptr_BR = log(ptr_sng * Mega_sum +1),
         mic_BR = log(mic_sng * Micro_sum +1),
         mol_BR = log(mol_sng * Molo_sum +1))

dat.2 <- dat.br %>% group_by(cell) %>%
  mutate(ptr_BR_1 = wrap(ptr_BR, n=1, order_by = month),
         mic_BR_1 = wrap(mic_BR, n=1, order_by = month),
         mol_BR_1 = wrap(mol_BR, n=1, order_by = month),
         ptr_BR_2 = wrap(ptr_BR, n=2, order_by = month),
         mic_BR_2 = wrap(mic_BR, n=2, order_by = month),
         mol_BR_2 = wrap(mol_BR, n=2, order_by = month),
         ptr_BR_3 = wrap(ptr_BR, n=3, order_by = month),
         mic_BR_3 = wrap(mic_BR, n=3, order_by = month),
         mol_BR_3 = wrap(mol_BR, n=3, order_by = month),
         ptr_BR_4 = wrap(ptr_BR, n=4, order_by = month),
         mic_BR_4 = wrap(mic_BR, n=4, order_by = month),
         mol_BR_4 = wrap(mol_BR, n=4, order_by = month),
         ptr_BR_5 = wrap(ptr_BR, n=5, order_by = month),
         mic_BR_5 = wrap(mic_BR, n=5, order_by = month),
         mol_BR_5 = wrap(mol_BR, n=5, order_by = month),
         ptr_BR_6 = wrap(ptr_BR, n=6, order_by = month),
         mic_BR_6 = wrap(mic_BR, n=6, order_by = month),
         mol_BR_6 = wrap(mol_BR, n=6, order_by = month), 
         OB_an_1 = wrap(OB_ann0_, n=1, order_by = month)) %>%
  ungroup %>%
  mutate(logPop = log(popDen + 1),
         NB_lDiv = log((Mam_sum - (Mega_sum + Molo_sum + Micro_sum))+1))
  ####DoubleMonthFrame####
dat.dbl <- dat %>%
  mutate(ptr_BR = log(ptr_dbl * Mega_sum +1),
         mic_BR = log(mic_dbl * Micro_sum +1),
         mol_BR = log(mol_dbl * Molo_sum +1))

dat.dbl <- dat.dbl %>% group_by(cell) %>%
  mutate(ptr_BR_1 = wrap(ptr_BR, n=1, order_by = month),
         mic_BR_1 = wrap(mic_BR, n=1, order_by = month),
         mol_BR_1 = wrap(mol_BR, n=1, order_by = month),
         ptr_BR_2 = wrap(ptr_BR, n=2, order_by = month),
         mic_BR_2 = wrap(mic_BR, n=2, order_by = month),
         mol_BR_2 = wrap(mol_BR, n=2, order_by = month),
         ptr_BR_3 = wrap(ptr_BR, n=3, order_by = month),
         mic_BR_3 = wrap(mic_BR, n=3, order_by = month),
         mol_BR_3 = wrap(mol_BR, n=3, order_by = month),
         ptr_BR_4 = wrap(ptr_BR, n=4, order_by = month),
         mic_BR_4 = wrap(mic_BR, n=4, order_by = month),
         mol_BR_4 = wrap(mol_BR, n=4, order_by = month),
         ptr_BR_5 = wrap(ptr_BR, n=5, order_by = month),
         mic_BR_5 = wrap(mic_BR, n=5, order_by = month),
         mol_BR_5 = wrap(mol_BR, n=5, order_by = month),
         ptr_BR_6 = wrap(ptr_BR, n=6, order_by = month),
         mic_BR_6 = wrap(mic_BR, n=6, order_by = month),
         mol_BR_6 = wrap(mol_BR, n=6, order_by = month), 
         OB_an_1 = wrap(OB_ann0_, n=1, order_by = month)) %>%
  ungroup %>%
  mutate(logPop = log(popDen + 1),
         NB_lDiv = log((Mam_sum - (Mega_sum + Molo_sum + Micro_sum))+1))

#### Create Owin object ####
library(raster) ; library(maptools) ; library(data.table)
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- raster(file.path(data.source, "CropMask.tif")) #used later as well
rf.poly <- rasterToPolygons(rf, fun = function(x){x==1}, dissolve = T)

library(spatstat)
regions <- slot(rf.poly, "polygons")
regions <- lapply(regions, function(x) { SpatialPolygons(list(x)) })
windows <- lapply(regions, as.owin)
#clean
rm(rf.poly, regions)


#### HUMAN SPILLOVER ####
  #### Single Month ####
    #### creating point object ####

humOB.W <- pppWeights(ob.col = OB_hum0_,
                      dataFrame = dat.2,
                      rasterGrid = rf)
    #### Covariates ####
hum.CoV <- c( "ptr_BR", "mic_BR", "mol_BR",
              "ptr_BR_1", "mic_BR_1", "mol_BR_1",
              "ptr_BR_2", "mic_BR_2", "mol_BR_2",
              "ptr_BR_3", "mic_BR_3", "mol_BR_3",
              "ptr_BR_4", "mic_BR_4", "mol_BR_4",
              "ptr_BR_5", "mic_BR_5", "mol_BR_5",
              "logPop", "OB_ann0_", "NB_lDiv",
              "fragIndex", "OB_an_1","month",
              "OB_hum0_",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
humW.df <- weightedDf(ob.col = OB_hum0_,
                      obWeights.df = humOB.W,
                      dataFrame = dat.2,
                      cols = hum.CoV)
humW.v <- humW.df$.mpl.W
    #### Model ####
hum.form <- as.formula(paste(".mpl.Y ","~", paste(hum.CoV[1:(length(hum.CoV)-5)],collapse = "+")))

hum.mod <- glm(hum.form,
               family=quasi(link="log", variance="mu"),
               weights=humW.v,
               data=humW.df)
summary(hum.mod)

    #### Predict ####

hum.pred <- predict.glm(hum.mod,newdata = humW.df, na.action = na.pass)
hum.pred.df <- cbind(hum.pred, humW.df)

out <- list()
for( i in 1:12){
  empty.raster <- rf
  empty.raster[] <- NA_real_
  cell.dance <- hum.pred.df %>%
    filter(month == i)
  empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$hum.pred
  names(empty.raster) <- paste0("pred_",i)
  out[[i]] <- empty.raster
}
lapply(out,plot)

  #### Double Month ###
  #### Double Month ####
    #### creating point object ####

humOB.W <- pppWeights(ob.col = OB_hum0_,
                      dataFrame = dat.dbl,
                      rasterGrid = rf)
    #### Covariates ####
hum.CoV <- c( "ptr_BR", "mic_BR", "mol_BR",
              "ptr_BR_1", "mic_BR_1", "mol_BR_1",
              "ptr_BR_3", "mic_BR_3", "mol_BR_3",
              "ptr_BR_5", "mic_BR_5", "mol_BR_5",
              "logPop", "OB_ann0_", "NB_lDiv",
              "fragIndex", "OB_an_1","month",
              "OB_hum0_",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
humW.df <- weightedDf(ob.col = OB_hum0_,
                      obWeights.df = humOB.W,
                      dataFrame = dat.dbl,
                      cols = hum.CoV)
humW.v <- humW.df$.mpl.W

    #### Model ####
hum.form <- as.formula(paste(".mpl.Y ","~", paste(hum.CoV[1:(length(hum.CoV)-5)],collapse = "+")))

hum.dbl <- glm(hum.form,
               family=quasi(link="log", variance="mu"),
               weights=humW.v,
               data=humW.df)
summary(hum.dbl)

    #### Predict ####

hum.pred.dbl <- predict.glm(hum.dbl,newdata = humW.df, na.action = na.pass)
hum.pred.df.dbl <- cbind(hum.pred.dbl, humW.df)

out.dbl <- list()
for( i in 1:12){
  empty.raster <- rf
  empty.raster[] <- NA_real_
  cell.dance <- hum.pred.df.dbl %>%
    filter(month == i)
  empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$hum.pred
  names(empty.raster) <- paste0("pred_",i)
  out.dbl[[i]] <- empty.raster
}
lapply(out.dbl,plot)
