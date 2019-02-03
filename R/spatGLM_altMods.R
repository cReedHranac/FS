#### spatGLM_altMods ####
## updated version of the spatGLM_qAIC scrip to include alternative models


source("R/helperFunctions.R")
library(data.table); library(tidyverse); library(rlang); library(lazyeval)
library(broom); library(gtools)
## Dave run this!
#clean.dir <- "."
#data.source <- "."

#### Functions ####
pppWeights <- function(ob.col, dataFrame, rasterGrid, owin = windows){
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

spatGLM <- function(ob.col, coV.v, dat, rGrid = rf,
                    project = F, save.name = NULL, leave.one.out = F){
  ###Function for applying the hybrid spatGLM to data
  ###Arguments
  ### ob.col <- outbreak column from dataframe
  ### coV.v <- vector of covariate names from dat
  # Note: last 5 items should include month, ob.col, x,y, and cell
  ### dat <- dataframe containing the relivent information
  ### rGrid <- rasterGrid from pppWeights obj (defalut to rf mask)
  ### project <- logical if the NoANimal projections should be made and written
  ### save.name <- name code to write files to 
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
    empty.raster[as.numeric(substring(cell.dance$cell,2))] <- exp(cell.dance$pred)
    names(empty.raster) <- paste0("pred_",i)
    rast[[i]] <- empty.raster
  }
  
  ### return w/o preojection
  items.out <- list(mod, pred.df, do.call(stack,rast))
  
  if(project == T){
    a <- spatGLM.AnimalMod(ob.col = !!ob.col, 
                           coV.v = coV.v,
                           rGrid = rGrid,
                           dat = dat)
    mod.stk <- a[[3]]
    out.loc <- file.path(mod.out.nov, save.name, "noAn")
    ifelse(!dir.exists(out.loc),
           dir.create(out.loc, recursive = T),
           F)
    writeRaster(mod.stk,
                out.loc, 
                format = "GTiff",
                bylayer = T,
                suffix = "numbers",
                overwrite= T)
    items.out <- list(model = mod,
                      pred.df = pred.df,
                      raster.projection = do.call(stack, rast),
                      Nonstocastic.projection = mod.stk)
  }
  
  if(leave.one.out == T){
    ## break appart to each individual case
    outbreak.obs <- which(mod$data$.mpl.Y > 0)
    W.df$Chunk <- NA
    ## maintain proportions of case to contorl
    W.df$Chunk[outbreak.obs] <- seq(1:32)
    things <- which(mod$data$.mpl.Y == 0)
    W.df$Chunk[things] <- rep_len(seq(1:32),length(things))
    
    ##begin
    roc.est <- list()
    for(i in 1:length(outbreak.obs)){
      j <- which(W.df$Chunk == i)
      ob <- outbreak.obs[[i]]
      
      y.hat <- predict.glm(update(mod, 
                                  data = mod$data[-j,],
                                  weights = W.v[-j]),
                           na.action = na.pass,
                           weights = W.v[-j])
      prob.0 <- gtools::inv.logit(y.hat)
      
      mod.dat <- mod$data[-j,][as.numeric(names(y.hat)),]
      roc.loo <- roc(mod.dat$OB_hum_imp, prob.0)
      roc.est[[i]] <- roc.loo$auc
      
    }
    mean.roc <- mean(unlist(roc.est))
    items.out <- c(items.out,
                   mean.roc = mean.roc)
  }
  
  return(items.out)
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

dfun <- function(object){
  ## Function from Ben Bolker to extract the overdispersion 
  with(object, sum((weights * residuals^2)[weights > 0])/df.residual)
}

qAIC <- function(x){
  ## Function for calculating the qAIC by hand since quasi models are funky
  # number of params
  k <- length(names(x[[1]]$coefficients))
  # overdispersion
  c.hat <- dfun(x[[1]])
  
  qAIC <- x[[1]]$deviance * c.hat + 2 * k
  
  out <- c(qAIC, c.hat, k)
  names(out) <- c("qAIC", "c.hat", "k")
  return(out)
}

resTabSimple <- function(x){
  ### FUnction for creating an easy results table 
  results <- cbind(tidy(x[[1]])[,-4], confint.default(x[[1]])[,1:2])
  tidy.table <- cbind(results[,1], round(results[,2:ncol(results)], 2))
  rownames(tidy.table) <- c()
  return(tidy.table)
}

spatGLM.AnimalMod <- function(ob.col, coV.v, dat, rGrid = rf){
  ###Function for applying the hybrid spatGLM to data
  ###Arguments
  ### ob.col <- outbreak column from dataframe
  ### coV.v <- vector of covariate names from dat
  # Note: last 5 items should include month, ob.col, x,y, and cell
  ### dat <- dataframe containing the relivent information
  ### rGrid <- rasterGrid from pppWeights obj (defalut to rf mask)
  
  ## Modified function to set animal outbreak beta to 0 for testing
  
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
  ## Modification ##
  W.df$OB_ann_imp_1 <- 0
  W.df$OB_ann_imp <- 0
  ##            ##
  pred <- predict.glm(mod,newdata = W.df, na.action = na.pass)
  pred.df <- cbind(pred, W.df)
  
  rast <- list()
  for( i in 1:12){
    empty.raster <- rGrid
    empty.raster[] <- NA_real_
    cell.dance <- pred.df %>%
      filter(month == i)
    empty.raster[as.numeric(substring(cell.dance$cell,2))] <- exp(cell.dance$pred)
    names(empty.raster) <- paste0("pred_",i)
    rast[[i]] <- empty.raster
  }
  
  items.out <- list(mod, pred.df, do.call(stack,rast))
  
}

#### Data ####
dat <- tbl_df(fread(file.path(mod.out.nov,"fullLongMaster.csv")))
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

#### Human Null Model (hNull) ####

h.null <- spatGLM(ob.col = OB_hum_imp,
                  coV.v = c("lPop", "lFrag", "hdl", "lnBm.div",
                            "OB_ann_imp", "OB_ann_imp_1",
                            "OB_hum_imp",
                            "month",  "x", "y", "cell"),
                  dat = dat,
                  project = T,
                  save.name = "h_null",
                  leave.one.out = T)
summary(h.null[[1]])
h.null.table <- resTabSimple(h.null)
write.csv(h.null.table, 
          "data/res_h_Null.csv",
          row.names = F)
#### Human Null Diversity Model (nDiv)####
h.nDiv <- spatGLM(ob.col = OB_hum_imp,
                  coV.v = c("lPop", "lFrag", "hdl", "lnBm.div",
                            "lBatDiv","OB_ann_imp", "OB_ann_imp_1",
                            "OB_hum_imp",
                            "month",  "x", "y", "cell"),
                  project = T,
                  save.name = "h_nDiv",
                  dat = dat,
                  leave.one.out = T)
summary(h.nDiv[[1]])
h.nDiv.table <- resTabSimple(h.nDiv)
write.csv(h.nDiv.table, 
          "data/res_h_nDiv.csv",
          row.names = F)

#### Human Null Split Bat Diversity Model nSBD (null Split Bat Diversity) ####
h.nSBD <- spatGLM(ob.col = OB_hum_imp,
                  coV.v = c("lPop", "lFrag", "hdl", "lnBm.div",
                            "lptrDiv", "lmicDiv", "lmolDiv",
                            "OB_ann_imp", "OB_ann_imp_1",
                            "OB_hum_imp",
                            "month",  "x", "y", "cell"),
                  project = T,
                  save.name = "h_nSBD",
                  dat = dat,
                  leave.one.out = T)
summary(h.nSBD[[1]])
h.nSBD.table <- resTabSimple(h.nSBD)
write.csv(h.nSBD.table, 
          "data/res_h_nSBD.csv",
          row.names = F)

#### Human Conditional Model A (No diversity) ####
h.cNDiv <- spatGLM(ob.col = OB_hum_imp,
                  coV.v = c("lBB.cond", "lBB.cond_2", "lBB.cond_4", "lBB.cond_6",
                            "lPop", "lFrag", "hdl", "lnBm.div",
                            "OB_ann_imp", "OB_ann_imp_1",
                            "OB_hum_imp",
                            "month",  "x", "y", "cell"),
                  project = T,
                  save.name = "h_cNDiv",
                  dat = dat,
                  leave.one.out = T)
summary(h.cNDiv[[1]])
h.cNDiv.table <- resTabSimple(h.cNDiv)
write.csv(h.cNDiv.table, 
          "data/res_h_cNDiv.csv",
          row.names = F)
#### Human Conditional Model B (with bat diversity) ####
h.cBDiv <- spatGLM(ob.col = OB_hum_imp,
                   coV.v = c("lBB.cond", "lBB.cond_2", "lBB.cond_4", "lBB.cond_6",
                             "lPop", "lFrag", "hdl", "lnBm.div", "lBatDiv",
                             "OB_ann_imp", "OB_ann_imp_1",
                             "OB_hum_imp",
                             "month",  "x", "y", "cell"),
                   project = T,
                   save.name = "h_cBDiv",
                   dat = dat,
                   leave.one.out = T)
summary(h.cBDiv[[1]])
h.cBDiv.table <- resTabSimple(h.cBDiv)
write.csv(h.cBDiv.table, 
          "data/res_h_cBDiv.csv",
          row.names = F)
#### Human Conditional Model B (Diversity Product) ####
h.cPDiv <- spatGLM(ob.col = OB_hum_imp,
                   coV.v = c("lBB.condDIV", "lBB.condDIV_2", "lBB.condDIV_4", "lBB.condDIV_6",
                             "lPop", "lFrag", "hdl", "lnBm.div",
                             "OB_ann_imp", "OB_ann_imp_1",
                             "OB_hum_imp",
                             "month",  "x", "y", "cell"),
                   project = T, 
                   save.name = "h_cPDiv",
                   dat = dat,
                   leave.one.out = T)
summary(h.cPDiv[[1]])
h.cPDiv.table <- resTabSimple(h.cPDiv)
write.csv(h.cPDiv.table, 
          "data/res_h_cPDiv.csv",
          row.names = F)
#### Human Force of Birthing ####
h.ORG <- spatGLM(ob.col = OB_hum_imp,
                   coV.v = c("lptr_BR", "lmic_BR", "lmol_BR",
                             "lptr_BR_2", "lmic_BR_2", "lmol_BR_2",
                             "lptr_BR_4", "lmic_BR_4", "lmol_BR_4",
                             "lptr_BR_6", "lmic_BR_6", "lmol_BR_6",
                             "lPop", "lFrag", "hdl", "lnBm.div",
                             "OB_ann_imp", "OB_ann_imp_1",
                             "OB_hum_imp",
                             "month",  "x", "y", "cell"), 
                 project = T,
                 save.name = "h_ORG",
                  dat = dat,
                 leave.one.out = T)
summary(h.ORG[[1]])
h.ORG.table <- resTabSimple(h.ORG)
write.csv(h.ORG.table, 
          "data/res_h_ORG.csv",
          row.names = F)
#### Human Breeding Probabilty + Diversity ####
h.Prb<- spatGLM(ob.col = OB_hum_imp,
                 coV.v = c("lptr_prob", "lmic_prob", "lmol_prob",
                           "lptr_prob_2", "lmic_prob_2", "lmol_prob_2",
                           "lptr_prob_4", "lmic_prob_4", "lmol_prob_4",
                           "lptr_prob_6", "lmic_prob_6", "lmol_prob_6",
                           "lptrDiv", "lmicDiv", "lmolDiv",
                           "lPop", "lFrag", "hdl", "lnBm.div",
                           "OB_ann_imp", "OB_ann_imp_1",
                           "OB_hum_imp",
                           "month",  "x", "y", "cell"),
                project = T,
                save.name = "h_Prb",
                dat = dat,
                leave.one.out = T)
summary(h.Prb[[1]])
h.Prb.table <- resTabSimple(h.Prb)
write.csv(h.Prb.table, 
          "data/res_h_Prb.csv",
          row.names = F)

#### HUman qAIC ####
h.mods <- list(h.null, h.nDiv, h.nSBD, h.cNDiv, h.cBDiv, h.cPDiv, h.ORG, h.Prb)
h.qAIC <- as.data.frame(do.call(rbind, lapply(h.mods, qAIC)),
                        row.names = c("null", "nDiv", "nSBD",
                                      "cNDiv", "cBDiv", "cPDiv", "ORg", "Prb"))
write.csv(h.qAIC, "data/humqAIC.csv")

#### ROC ####
r <- list()
for(i in 1:length(h.mods)){
  r[[i]] <- mean(h.mods[[i]]$mean.roc)
}
roc.df <- as.data.frame(do.call(rbind, r),
                        row.names = c("null", "nDiv", "nSBD",
                                      "cNDiv", "cBDiv", "cPDiv", "ORg", "Prb"))
write.csv(roc.df, "data/humROC.csv")
#### Cross Validation #### 
# spatGLM.loo <- function(ob.col, coV.v, dat, rGrid = rf){
#   ###Function for applying the hybrid spatGLM to data
#   ###Arguments
#   ### ob.col <- outbreak column from dataframe
#   ### coV.v <- vector of covariate names from dat
#   # Note: last 5 items should include month, ob.col, x,y, and cell
#   ### dat <- dataframe containing the relivent information
#   ### rGrid <- rasterGrid from pppWeights obj (defalut to rf mask)
#   #### creating point object ####
#   ob.col <- enquo(ob.col)
#   ppp.ob <- pppWeights(ob.col = UQ(ob.col),
#                        dataFrame = dat,
#                        rasterGrid = rGrid)
#   #### Weighted dataframe/ weights vector ####
#   W.df <- weightedDf(ob.col = UQ(ob.col),
#                      obWeights.df = ppp.ob,
#                      dataFrame = dat,
#                      cols = coV.v)
#   W.v <- W.df$.mpl.W #weights vector
#   #### Model ####
#   form <- as.formula(paste(".mpl.Y ","~",
#                            paste(coV.v[1:(length(coV.v)-5)],collapse = "+")))
#   
#   mod <- glm(form,
#              family=quasi(link="log", variance="mu"),
#              weights=W.v,
#              data=W.df)
#   
#   #### Leave one out ####
#   ## Atempt 1: 
#   ## asign each case a # between 1:32, then split the data into 32nds to retain the
#   ## portion of cases to non-cases
#   
#   outbreak.obs <- which(mod$data$.mpl.Y > 0)
#   W.df$Chunk <- NA
#   W.df$Chunk[outbreak.obs] <- seq(1:32)
#   things <- which(mod$data$.mpl.Y == 0)
#   W.df$Chunk[things] <- rep_len(seq(1:32),length(things))
#   
#   
#   
#   roc.loo <- list()
#   roc.est <- list()
#   for(i in 1:length(outbreak.obs)){
#     j <- which(W.df$Chunk == i)
#     ob <- outbreak.obs[[i]]
#     
#     y.hat <- predict.glm(update(mod, 
#                                 data = mod$data[-j,],
#                                  weights = W.v[-j]),
#                           na.action = na.pass,
#                           weights = W.v[-j])
#     prob.0 <- gtools::inv.logit(y.hat)
#     
#     mod.dat <- mod$data[-j,][as.numeric(names(y.hat)),]
#     roc.loo[[i]] <- roc(mod.dat$OB_hum_imp, prob.0)
#     roc.est[[i]] <- roc.loo[[i]]$auc
#     
#   }
#   roc.scores <- do.call(rbind, roc.est)
#   
#   # items.out <- list(mod, pred.df, rast)
#   items.out <- list(mod, roc.scores)
# }
# 
# h.nSBD.loo <- spatGLM.loo(ob.col = OB_hum_imp,
#                   coV.v = c("lPop", "lFrag", "hdl", "lnBm.div",
#                             "lptrDiv", "lmicDiv", "lmolDiv",
#                             "OB_ann_imp", "OB_ann_imp_1",
#                             "OB_hum_imp",
#                             "month",  "x", "y", "cell"),
#                   dat = dat)
# h.Prb.loo <- spatGLM.loo(ob.col = OB_hum_imp,
#                 coV.v = c("lptr_prob", "lmic_prob", "lmol_prob",
#                           "lptr_prob_2", "lmic_prob_2", "lmol_prob_2",
#                           "lptr_prob_4", "lmic_prob_4", "lmol_prob_4",
#                           "lptr_prob_6", "lmic_prob_6", "lmol_prob_6",
#                           "lptrDiv", "lmicDiv", "lmolDiv",
#                           "lPop", "lFrag", "hdl", "lnBm.div",
#                           "OB_ann_imp", "OB_ann_imp_1",
#                           "OB_hum_imp",
#                           "month",  "x", "y", "cell"),
#                 dat = dat)
# mean(h.nSBD.loo[[2]]) 0.5505365
# mean(h.Prb.loo[[2]]) 0.5447542
