#### spatGLM clean with imputed data ####
source("R/helperFunctions.R")
library(data.table); library(dplyr); library(tidyr); library(rlang); library(lazyeval)

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
  obGrid <- obWeights.df %>% ##From the raster grid
    mutate_(.dots = setNames(list(interp(~gridMe(x, mx2),
                                 x=as.name("x"),
                                 mx2 = min(dataFrame[,"x"])),
                                 interp(~gridMe(y, my2),
                                        y = as.name("y"),
                                        my2 = min(dataFrame[,"y"]))),
                                nm = c("x_grid", "y_grid")))
      
  datGrid <- dat.sub %>% ##From the ppp grid
    mutate_(.dots = setNames(list(interp(~gridMe(x, mx2),
                                         x=as.name("x"),
                                         mx2 = min(obWeights.df[,"x"])),
                                  interp(~gridMe(y, my2),
                                         y = as.name("y"),
                                         my2 = min(obWeights.df[,"y"]))),
                                  nm = c("x_grid", "y_grid")))%>%
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

hum.sng.raw.spatGLM <- spatGLM(ob.col = OB_hum_imp,
                                  coV.v = c( "ptr_sng_raw_BR", "mic_sng_raw_BR", "mol_sng_raw_BR",
                                             "ptr_sng_raw_BR_1", "mic_sng_raw_BR_1", "mol_sng_raw_BR_1",
                                             "ptr_sng_raw_BR_2", "mic_sng_raw_BR_2", "mol_sng_raw_BR_2",
                                             "ptr_sng_raw_BR_3", "mic_sng_raw_BR_3", "mol_sng_raw_BR_3",
                                             "ptr_sng_raw_BR_4", "mic_sng_raw_BR_4", "mol_sng_raw_BR_4",
                                             "ptr_sng_raw_BR_5", "mic_sng_raw_BR_5", "mol_sng_raw_BR_5",
                                             "logPop", "OB_ann_imp", "NB_lDiv",
                                             "fragIndex", "OB_ann_imp_1","month",
                                             "OB_hum_imp",  "x", "y", "cell"),
                                  dat = dat)
#### No Micro ####
hum.sng.raw.NM.spatGLM <- spatGLM(ob.col = OB_hum_imp,
                               coV.v = c( "ptr_sng_raw_BR", "mol_sng_raw_BR",
                                          "ptr_sng_raw_BR_1", "mol_sng_raw_BR_1",
                                          "ptr_sng_raw_BR_2", "mol_sng_raw_BR_2",
                                          "ptr_sng_raw_BR_3", "mol_sng_raw_BR_3",
                                          "ptr_sng_raw_BR_4", "mol_sng_raw_BR_4",
                                          "ptr_sng_raw_BR_5", "mol_sng_raw_BR_5",
                                          "logPop", "OB_ann_imp", "NB_lDiv",
                                          "fragIndex", "OB_ann_imp_1","month",
                                          "OB_hum_imp",  "x", "y", "cell"),
                               dat = dat)
summary(hum.sng.raw.NM.spatGLM[[1]])

#### Imp ####
hum.sng.imp.spatGLM <- spatGLM(ob.col =  OB_hum_imp, 
                               coV.v = c( "ptr_sng_imp_BR", "mic_sng_imp_BR", "mol_sng_imp_BR",
                                          "ptr_sng_imp_BR_1", "mic_sng_imp_BR_1", "mol_sng_imp_BR_1",
                                          "ptr_sng_imp_BR_2", "mic_sng_imp_BR_2", "mol_sng_imp_BR_2",
                                          "ptr_sng_imp_BR_3", "mic_sng_imp_BR_3", "mol_sng_imp_BR_3",
                                          "ptr_sng_imp_BR_4", "mic_sng_imp_BR_4", "mol_sng_imp_BR_4",
                                          "ptr_sng_imp_BR_5", "mic_sng_imp_BR_5", "mol_sng_imp_BR_5",
                                          "logPop", "OB_ann_imp", "NB_lDiv",
                                          "fragIndex", "OB_ann_imp_1","month",
                                          "OB_hum_imp",  "x", "y", "cell"),
                               dat = dat)
summary(hum.sng.imp.spatGLM[[1]])

#### No Micro ####
hum.sng.imp.NM.spatGLM <- spatGLM(ob.col =  OB_hum_imp, 
                               coV.v = c( "ptr_sng_imp_BR", "mol_sng_imp_BR",
                                          "ptr_sng_imp_BR_1", "mol_sng_imp_BR_1",
                                          "ptr_sng_imp_BR_2", "mol_sng_imp_BR_2",
                                          "ptr_sng_imp_BR_3", "mol_sng_imp_BR_3",
                                          "ptr_sng_imp_BR_4", "mol_sng_imp_BR_4",
                                          "ptr_sng_imp_BR_5", "mol_sng_imp_BR_5",
                                          "logPop", "OB_ann_imp", "NB_lDiv",
                                          "fragIndex", "OB_ann_imp_1","month",
                                          "OB_hum_imp",  "x", "y", "cell"),
                               dat = dat)
summary(hum.sng.imp.NM.spatGLM[[1]])

  #### Double Month ####
#### Raw #### 
    #### creating point object ####
hum.dbl.raw.spatGLM <- spatGLM(OB_hum_imp,
                               coV.v = c( "ptr_dbl_raw_BR", "mic_dbl_raw_BR", "mol_dbl_raw_BR",
                                          "ptr_dbl_raw_BR_2", "mic_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
                                          "ptr_dbl_raw_BR_4", "mic_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
                                          "ptr_dbl_raw_BR_6", "mic_dbl_raw_BR_6", "mol_dbl_raw_BR_6",
                                          "logPop", "OB_ann_imp", "NB_lDiv",
                                          "fragIndex", "OB_ann_imp_1","month",
                                          "OB_hum_imp",  "x", "y", "cell"),
                               dat = dat)
summary(hum.dbl.raw.spatGLM[[1]])

#### No Micro
hum.dbl.raw.NM.spatGLM <- spatGLM(OB_hum_imp,
                               coV.v = c( "ptr_dbl_raw_BR", "mol_dbl_raw_BR",
                                          "ptr_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
                                          "ptr_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
                                          "ptr_dbl_raw_BR_6",  "mol_dbl_raw_BR_6",
                                          "logPop", "OB_ann_imp", "NB_lDiv",
                                          "fragIndex", "OB_ann_imp_1","month",
                                          "OB_hum_imp",  "x", "y", "cell"),
                               dat = dat)
summary(hum.dbl.raw.NM.spatGLM[[1]])


#### Imp #### 
hum.dbl.imp.spatGLM <-spatGLM(ob.col = OB_hum_imp,
                              coV.v = c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
                                         "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                         "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                         "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                         "logPop", "OB_ann_imp", "NB_lDiv",
                                         "fragIndex", "OB_ann_imp_1","month",
                                         "OB_hum_imp",  "x", "y", "cell"),
                              dat= dat)

summary(hum.dbl.imp.spatGLM[[1]])
#### No Micro
hum.dbl.imp.NM.spatGLM <-spatGLM(ob.col = OB_hum_imp,
                              coV.v = c( "ptr_dbl_imp_BR",  "mol_dbl_imp_BR",
                                         "ptr_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                         "ptr_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                         "ptr_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                         "logPop", "OB_ann_imp", "NB_lDiv",
                                         "fragIndex", "OB_ann_imp_1","month",
                                         "OB_hum_imp",  "x", "y", "cell"),
                              dat= dat)

summary(hum.dbl.imp.NM.spatGLM[[1]])


#### Animal Outbreaks ####
  #### Double Month ####
#### Imp #### 
# it does have imputed animal, and human outbreak occurence points though. it doesn't work without

ann.dbl.imp.spatGLM <- spatGLM(ob.col = OB_ann_imp,
                               coV.v = c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
                                          "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                          "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                          "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                          "NB_lDiv",
                                          "fragIndex", "month",
                                          "OB_ann_imp",  "x", "y", "cell"),
                               dat = dat)
summary(ann.dbl.imp.spatGLM[[1]])

#### No Micro
ann.dbl.imp.NM.spatGLM <- spatGLM(ob.col = OB_ann_imp,
                               coV.v = c( "ptr_dbl_imp_BR",  "mol_dbl_imp_BR",
                                          "ptr_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                          "ptr_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                          "ptr_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                          "NB_lDiv",
                                          "fragIndex", "month",
                                          "OB_ann_imp",  "x", "y", "cell"),
                               dat = dat)
summary(ann.dbl.imp.NM.spatGLM[[1]])

#### Raw ####
ann.dbl.raw.spatGLM <- spatGLM(ob.col = OB_ann_imp,
                               coV.v =  c( "ptr_dbl_raw_BR", "mic_dbl_raw_BR", "mol_dbl_raw_BR",
                                           "ptr_dbl_raw_BR_2", "mic_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
                                           "ptr_dbl_raw_BR_4", "mic_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
                                           "ptr_dbl_raw_BR_6", "mic_dbl_raw_BR_6", "mol_dbl_raw_BR_6",
                                           "NB_lDiv",
                                           "fragIndex", "month",
                                           "OB_ann_imp",  "x", "y", "cell"),
                               dat = dat)
summary(ann.dbl.raw.spatGLM[[1]])

#### No Micro
ann.dbl.raw.NM.spatGLM <- spatGLM(ob.col = OB_ann_imp,
                               coV.v =  c( "ptr_dbl_raw_BR", "mol_dbl_raw_BR",
                                           "ptr_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
                                           "ptr_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
                                           "ptr_dbl_raw_BR_6", "mol_dbl_raw_BR_6",
                                           "NB_lDiv",
                                           "fragIndex", "month",
                                           "OB_ann_imp",  "x", "y", "cell"),
                               dat = dat)
summary(ann.dbl.raw.NM.spatGLM[[1]])
