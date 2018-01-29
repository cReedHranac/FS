#### spatGLM clean with imputed data ####
source("R/helperFunctions.R")
library(data.table); library(dplyr); library(tidyr)

## Dave run this!
#clean.dir <- "."
#data.source <- "."

#### Functions ####
pppWeights <- function(ob.col, dataFrame, rasterGrid, owin = windows){
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
            window = owin[[1]])
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
    mutate_(x_grid := ~gridMe(x, min(dataFrame[,"x"])), ##From the raster grid
            y_grid := ~gridMe(y, min(dataFrame[,"y"])))
  datGrid <- dat.sub %>% 
    mutate_(x_grid := ~gridMe(x, min(obWeights.df[,"x"])), ##From the ppp grid
            y_grid := ~gridMe(y, min(obWeights.df[,"y"]))) %>%
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

#### Data ####
dat <- tbl_df(fread(file.path(clean.dir, "longMaster.csv")))
# library(skimr)
#  skim(dat)
dat <- dat %>%
  dplyr::group_by(cell) %>%
  dplyr::mutate(OB_ann_imp_1 = wrap(OB_ann_imp, 1, month)) %>%
  ungroup
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
#### Raw ####
    #### creating point object ####

hum.sng.raw.OB.W <- pppWeights(ob.col = OB_hum_imp,
                      dataFrame = dat,
                      rasterGrid = rf)
    #### Covariates ####
hum.sng.raw.CoV <- c( "ptr_sng_raw_BR", "mic_sng_raw_BR", "mol_sng_raw_BR",
              "ptr_sng_raw_BR_1", "mic_sng_raw_BR_1", "mol_sng_raw_BR_1",
              "ptr_sng_raw_BR_2", "mic_sng_raw_BR_2", "mol_sng_raw_BR_2",
              "ptr_sng_raw_BR_3", "mic_sng_raw_BR_3", "mol_sng_raw_BR_3",
              "ptr_sng_raw_BR_4", "mic_sng_raw_BR_4", "mol_sng_raw_BR_4",
              "ptr_sng_raw_BR_5", "mic_sng_raw_BR_5", "mol_sng_raw_BR_5",
              "logPop", "OB_ann_imp", "NB_lDiv",
              "fragIndex", "OB_ann_imp_1","month",
              "OB_hum_imp",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
humW.sng.raw.df <- weightedDf(ob.col = OB_hum_imp,
                      obWeights.df = hum.sng.raw.OB.W,
                      dataFrame = dat,
                      cols = hum.sng.raw.CoV)
humW.sng.raw.v <- humW.sng.raw.df$.mpl.W
    #### Model ####
hum.sng.raw.form <- as.formula(paste(".mpl.Y ","~",
                                     paste(hum.sng.raw.CoV[1:(length(hum.sng.raw.CoV)-5)],collapse = "+")))

hum.sng.raw.mod <- glm(hum.sng.raw.form,
               family=quasi(link="log", variance="mu"),
               weights=humW.sng.raw.v,
               data=humW.sng.raw.df)
summary(hum.sng.raw.mod)

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
#### Imp ####
    #### creating point object ####

humOB.sng.imp.W <- pppWeights(ob.col = OB_hum_imp,
                      dataFrame = dat,
                      rasterGrid = rf)
    #### Covariates ####
hum.sng.imp.CoV <- c( "ptr_sng_imp_BR", "mic_sng_imp_BR", "mol_sng_imp_BR",
              "ptr_sng_imp_BR_1", "mic_sng_imp_BR_1", "mol_sng_imp_BR_1",
              "ptr_sng_imp_BR_2", "mic_sng_imp_BR_2", "mol_sng_imp_BR_2",
              "ptr_sng_imp_BR_3", "mic_sng_imp_BR_3", "mol_sng_imp_BR_3",
              "ptr_sng_imp_BR_4", "mic_sng_imp_BR_4", "mol_sng_imp_BR_4",
              "ptr_sng_imp_BR_5", "mic_sng_imp_BR_5", "mol_sng_imp_BR_5",
              "logPop", "OB_ann_imp", "NB_lDiv",
              "fragIndex", "OB_ann_imp_1","month",
              "OB_hum_imp",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
humW.sng.imp.df <- weightedDf(ob.col = OB_hum_imp,
                      obWeights.df = humOB.sng.imp.W,
                      dataFrame = dat,
                      cols = hum.sng.imp.CoV)
humW.sng.imp.v <- humW.sng.imp.df$.mpl.W

    #### Model ####
hum.sng.imp.form <- as.formula(paste(".mpl.Y ","~", 
                                     paste(hum.sng.imp.CoV[1:(length(hum.sng.imp.CoV)-5)],collapse = "+")))

hum.sng.imp.mod <- glm(hum.sng.imp.form,
               family=quasi(link="log", variance="mu"),
               weights=humW.sng.imp.v,
               data=humW.sng.imp.df)
summary(hum.sng.imp.mod)

    #### Predict ####

hum.sng.imp.pred <- predict.glm(hum.sng.imp.mod,newdata = humW.sng.imp.df, na.action = na.pass)
hum.sng.imp.pred.df <- cbind(hum.sng.imp.pred, humW.sng.imp.df)

out.hum.sng.imp <- list()
for( i in 1:12){
  empty.raster <- rf
  empty.raster[] <- NA_real_
  cell.dance <- hum.sng.imp.pred.df %>%
    filter(month == i)
  empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$hum.sng.imp.pred
  names(empty.raster) <- paste0("pred_",i)
  out.hum.sng.imp[[i]] <- empty.raster
}
lapply(out.hum.sng.imp,plot)



  #### Double Month ####
#### Raw #### 
# it does have imputed animal, and human outbreak occurence points though. it doesn't work without
    #### creating point object ####

hum.dbl.raw.OB.W <- pppWeights(ob.col = OB_hum_imp,
                      dataFrame = dat,
                      rasterGrid = rf)
    #### Covariates ####
hum.dbl.raw.CoV <- c( "ptr_dbl_raw_BR", "mic_dbl_raw_BR", "mol_dbl_raw_BR",
              "ptr_dbl_raw_BR_2", "mic_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
              "ptr_dbl_raw_BR_4", "mic_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
              "ptr_dbl_raw_BR_6", "mic_dbl_raw_BR_6", "mol_dbl_raw_BR_6",
              "logPop", "OB_ann_imp", "NB_lDiv",
              "fragIndex", "OB_ann_imp_1","month",
              "OB_hum_imp",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
humW.dbl.raw.df <- weightedDf(ob.col = OB_hum_imp,
                      obWeights.df = hum.dbl.raw.OB.W,
                      dataFrame = dat,
                      cols = hum.dbl.raw.CoV)
humW.dbl.raw.v <- humW.dbl.raw.df$.mpl.W

    #### Model ####
hum.dbl.raw.form <- as.formula(paste(".mpl.Y ","~",
                                     paste(hum.dbl.raw.CoV[1:(length(hum.dbl.raw.CoV)-5)],collapse = "+")))

hum.dbl.raw.mod <- glm(hum.dbl.raw.form,
               family=quasi(link="log", variance="mu"),
               weights=humW.dbl.raw.v,
               data=humW.dbl.raw.df)
summary(hum.dbl.raw.mod)

    #### Predict ####

hum.pred.dbl.raw <- predict.glm(hum.dbl.raw.mod,newdata = humW.dbl.raw.df, na.action = na.pass)
hum.pred.df.dbl.raw <- cbind(hum.pred.dbl.raw, humW.dbl.raw.df)

out.hum.dbl.raw <- list()
for( i in 1:12){
  empty.raster <- rf
  empty.raster[] <- NA_real_
  cell.dance <- hum.pred.df.dbl.raw %>%
    filter(month == i)
  empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$hum.pred.dbl.raw
  names(empty.raster) <- paste0("pred_",i)
  out.hum.dbl.raw[[i]] <- empty.raster
}
lapply(out.hum.dbl.raw,plot)

#### Imp #### 
# it does have imputed animal, and human outbreak occurence points though. it doesn't work without
    #### creating point object ####

humOB.W <- pppWeights(ob.col = OB_hum_imp,
                      dataFrame = dat,
                      rasterGrid = rf)
    #### Covariates ####
hum.CoV <- c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
              "ptr_dbl_imp_BR_1", "mic_dbl_imp_BR_1", "mol_dbl_imp_BR_1",
              "ptr_dbl_imp_BR_3", "mic_dbl_imp_BR_3", "mol_dbl_imp_BR_3",
              "ptr_dbl_imp_BR_5", "mic_dbl_imp_BR_5", "mol_dbl_imp_BR_5",
              "logPop", "OB_ann_imp", "NB_lDiv",
              "fragIndex", "OB_ann_imp_1","month",
              "OB_hum_imp",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
humW.df <- weightedDf(ob.col = OB_hum_imp,
                      obWeights.df = humOB.W,
                      dataFrame = dat,
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




#### Animal Outbreaks ####
  #### Double Month ####
#### Imp #### 
# it does have imputed animal, and human outbreak occurence points though. it doesn't work without
    #### creating point object ####

annOB.dbl.imp.W <- pppWeights(ob.col = OB_ann_imp,
                      dataFrame = dat,
                      rasterGrid = rf)
    #### Covariates ####
ann.dbl.imp.CoV <- c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
              "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
              "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
              "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
              "logPop", "NB_lDiv",
              "fragIndex", "month",
              "OB_ann_imp",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
annW.dbl.imp.df <- weightedDf(ob.col = OB_ann_imp,
                      obWeights.df = annOB.dbl.imp.W,
                      dataFrame = dat,
                      cols = ann.dbl.imp.CoV)
annW.dbl.imp.v <- annW.dbl.imp.df$.mpl.W

    #### Model ####
ann.dbl.imp.form <- as.formula(paste(".mpl.Y ","~",
                                     paste(ann.dbl.imp.CoV[1:(length(ann.dbl.imp.CoV)-5)],collapse = "+")))

ann.dbl.imp.mod <- glm(ann.dbl.imp.form,
               family=quasi(link="log", variance="mu"),
               weights=annW.dbl.imp.v,
               data=annW.dbl.imp.df)
summary(ann.dbl.imp.mod)

    #### Predict ####

ann.pred.dbl.imp <- predict.glm(ann.dbl.imp.mod,newdata = annW.dbl.imp.df, na.action = na.pass)
ann.pred.df.dbl.imp <- cbind(ann.pred.dbl.imp, annW.dbl.imp.df)

out.ann.dbl.imp <- list()
for( i in 1:12){
  empty.raster <- rf
  empty.raster[] <- NA_real_
  cell.dance <- ann.pred.df.dbl.imp %>%
    filter(month == i)
  empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$ann.pred.dbl.imp
  names(empty.raster) <- paste0("pred_",i)
  out.ann.dbl.imp[[i]] <- empty.raster
}
lapply(out.ann.dbl.imp,plot)

#### Raw ####
    #### creating point object ####
annOB.dbl.raw.W <- pppWeights(ob.col = OB_ann_imp,
                      dataFrame = dat,
                      rasterGrid = rf)
    #### Covariates ####
ann.dbl.raw.CoV <- c( "ptr_dbl_raw_BR", "mic_dbl_raw_BR", "mol_dbl_raw_BR",
              "ptr_dbl_raw_BR_2", "mic_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
              "ptr_dbl_raw_BR_4", "mic_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
              "ptr_dbl_raw_BR_6", "mic_dbl_raw_BR_6", "mol_dbl_raw_BR_6",
              "logPop", "NB_lDiv",
              "fragIndex", "month",
              "OB_ann_imp",  "x", "y", "cell")

    #### Weighted dataframe/ weights vector ####
annW.dbl.raw.df <- weightedDf(ob.col = OB_ann_imp,
                      obWeights.df = annOB.dbl.raw.W,
                      dataFrame = dat,
                      cols = ann.dbl.raw.CoV)
annW.dbl.raw.v <- annW.dbl.raw.df$.mpl.W

    #### Model ####
ann.dbl.raw.form <- as.formula(paste(".mpl.Y ","~",
                                     paste(ann.dbl.raw.CoV[1:(length(ann.dbl.raw.CoV)-5)],collapse = "+")))

ann.dbl.raw.mod <- glm(ann.dbl.raw.form,
               family=quasi(link="log", variance="mu"),
               weights=annW.dbl.raw.v,
               data=annW.dbl.raw.df)
summary(ann.dbl.raw.mod)

    #### Predict ####

ann.pred.dbl.raw <- predict.glm(ann.dbl.raw.mod, newdata = annW.dbl.raw.df, na.action = na.pass)
ann.pred.df.dbl.raw <- cbind(ann.pred.dbl.raw, annW.dbl.raw.df)

out.ann.dbl.raw <- list()
for( i in 1:12){
  empty.raster <- rf
  empty.raster[] <- NA_real_
  cell.dance <- ann.pred.df.dbl.raw %>%
    filter(month == i)
  empty.raster[as.numeric(substring(cell.dance$cell,2))] <- cell.dance$ann.pred.dbl.raw
  names(empty.raster) <- paste0("pred_",i)
  out.ann.dbl.raw[[i]] <- empty.raster
}
lapply(out.ann.dbl.raw,plot)
