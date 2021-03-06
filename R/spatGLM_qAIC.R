#### spatGLM_qAIC ####

source("R/helperFunctions.R")
library(data.table); library(dplyr); library(tidyr); library(rlang); library(lazyeval)
library(broom)
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
    empty.raster[as.numeric(substring(cell.dance$cell,2))] <- exp(cell.dance$pred)
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
  
  items.out <- list(mod, pred.df, rast)
}

#### Data ####
dat <- tbl_df(fread(file.path(clean.dir.nov, "longMaster.csv")))
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

#### Human Model ####

#############################
### Full
hum.full <-spatGLM(ob.col = OB_hum_imp,
                              coV.v = c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
                                         "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                         "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                         "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                         "OB_ann_imp","OB_ann_imp_1",
                                         "hdl","logPop","lnBm.div","lFrag",
                                         "OB_hum_imp","month",
                                         "x", "y", "cell"),
                              dat= dat)

summary(hum.full[[1]])

# tidy and write out
humFullTable <- resTabSimple(hum.full)
write.csv(humFullTable, "data/HumSpGLMRes_Org.csv", row.names = F)

mod.stk <- do.call(stack, hum.full[[3]])
writeRaster(mod.stk, file.path(mod.out.nov, "SpGLMRes_Org", "hum"),format = "GTiff",
            bylayer = T, suffix = "numbers", overwrite = T)

## Remove animal bits for the outbreak fitting
hum.full.noAn <-spatGLM.AnimalMod(ob.col = OB_hum_imp,
                   coV.v = c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
                              "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                              "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                              "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                              "OB_ann_imp","OB_ann_imp_1",
                              "hdl","logPop","lnBm.div","lFrag",
                              "OB_hum_imp","month",
                              "x", "y", "cell"),
                   dat= dat)
mod.stk <- do.call(stack, hum.full.noAn[[3]])
writeRaster(mod.stk,
            file.path(mod.out.nov, "SpGLMRes_Org", "humNoAn"),
            format = "GTiff",
            bylayer=T,
            suffix = "numbers",
            overwrite = T)
##############################
## Null model
hum.null <-spatGLM(ob.col = OB_hum_imp,
                      coV.v = c(  "OB_ann_imp", "OB_ann_imp_1","lnBm.div","logPop",
                                 "lFrag", "month",
                                 "OB_hum_imp",  "x", "y", "cell"),
                      dat= dat)
summary(hum.null[[1]])

humNullTable <- resTabSimple(hum.null)
write.csv(humNullTable, "data/HumNullSpGLMRes_Orgcsv", row.names = F)


########################################
## Alternative mod (raw probablilities rather then multiplied with diversity)
hum.prob <- spatGLM(ob.col = OB_hum_imp,
                   coV.v = c( "ptr_dbl_imp", "mic_dbl_imp", "mol_dbl_imp",
                              "ptr_dbl_imp_Prob_2", "mic_dbl_imp_Prob_2", "mol_dbl_imp_Prob_2",
                              "ptr_dbl_imp_Prob_4", "mic_dbl_imp_Prob_4", "mol_dbl_imp_Prob_4",
                              "ptr_dbl_imp_Prob_6", "mic_dbl_imp_Prob_6", "mol_dbl_imp_Prob_6",
                              "OB_ann_imp","OB_ann_imp_1", "BatDiv",
                              "hdl","logPop","lnBm.div","lFrag",
                              "OB_hum_imp","month",
                              "x", "y", "cell"),
                   dat= dat)
summary(hum.prob[[1]])

# tidy and write out
humprobTable <- resTabSimple(hum.prob)
write.csv(humprobTable, "data/HumSpGLMRes_Prob.csv", row.names = F)

mod.stk <- do.call(stack, hum.prob[[3]])
writeRaster(mod.stk, file.path(mod.out.nov, "SpGLMRes_Prob", "hum"),format = "GTiff",
            bylayer = T, suffix = "numbers", overwrite = T)
## Modified to remove all animal outbreak information post model fit
hum.NoAn <-spatGLM.AnimalMod(ob.col = OB_hum_imp,
                             coV.v = c( "ptr_dbl_imp", "mic_dbl_imp", "mol_dbl_imp",
                                        "ptr_dbl_imp_Prob_2", "mic_dbl_imp_Prob_2", "mol_dbl_imp_Prob_2",
                                        "ptr_dbl_imp_Prob_4", "mic_dbl_imp_Prob_4", "mol_dbl_imp_Prob_4",
                                        "ptr_dbl_imp_Prob_6", "mic_dbl_imp_Prob_6", "mol_dbl_imp_Prob_6",
                                        "logPop", "OB_ann_imp", "lnBm.div", "BatDiv",
                                        "hdl","lFrag", "OB_ann_imp_1","month",
                                        "OB_hum_imp",  "x", "y", "cell"),
                             dat= dat)

summary(hum.NoAn[[1]])
humNoAnTable <- resTabSimple(hum.NoAn)
write.csv(humNoAnTable, "data/HumSpGLMRes_NoAn_Prob.csv", row.names = F)
modNoAn.stk <- do.call(stack, hum.NoAn[[3]])
writeRaster(modNoAn.stk,
            file.path(mod.out.nov, "SpGLMRes_Prob", "humNoAn"),
            format = "GTiff",
            bylayer = T,
            suffix = "numbers",
            overwrite = T)


###################################
## Alternative mod 2: totals
hum.sum <- spatGLM(ob.col = OB_hum_imp,
                    coV.v = c( "BB.cond",
                               "BB.cond_2",
                               "BB.cond_4",
                               "BB.cond_6",
                               "OB_ann_imp","OB_ann_imp_1", "BatDiv",
                               "hdl","logPop","lnBm.div","lFrag",
                               "OB_hum_imp","month",
                               "x", "y", "cell"),
                    dat= dat)
summary(hum.sum[[1]])

# tidy and write out
humsumTable <- resTabSimple(hum.sum)
write.csv(humsumTable, "data/HumSpGLMRes_sum.csv", row.names = F)

mod.stk <- do.call(stack, hum.sum[[3]])
writeRaster(mod.stk, file.path(mod.out.nov, "SpGLMRes_sum", "hum"),format = "GTiff",
            bylayer = T, suffix = "numbers", overwrite = T)

## No animal version 
hum.sumNoAn <- spatGLM.AnimalMod(ob.col = OB_hum_imp,
                   coV.v = c( "BB.cond",
                              "BB.cond_2",
                              "BB.cond_4",
                              "BB.cond_6",
                              "OB_ann_imp","OB_ann_imp_1", "BatDiv",
                              "hdl","logPop","lnBm.div","lFrag",
                              "OB_hum_imp","month",
                              "x", "y", "cell"),
                   dat= dat)

sumNoAn <- do.call(stack, hum.sumNoAn[[3]])
writeRaster(sumNoAn,
            file.path(mod.out.nov, "SpGLMRes_sum", "humNoAn"),
            format = "GTiff",
            bylayer = T,
            suffix = "numbers",
            overwrite = T)

####################################
## Alternative mod 3: totals x BatDiv
## modify dataframe
dat.mod <- dat %>%
  mutate(BB.condDIV = BB.cond * BatDiv,
         BB.condDIV_2 = BB.cond_2 * BatDiv,
         BB.condDIV_4 = BB.cond_4 * BatDiv,
         BB.condDIV_6 = BB.cond_6 * BatDiv)

hum.sumDIV <- spatGLM(ob.col = OB_hum_imp,
                   coV.v = c( "BB.condDIV",
                              "BB.condDIV_2",
                              "BB.condDIV_4",
                              "BB.condDIV_6",
                              "OB_ann_imp","OB_ann_imp_1", 
                              "hdl","logPop","lnBm.div","lFrag",
                              "OB_hum_imp","month",
                              "x", "y", "cell"),
                   dat= dat.mod)
summary(hum.sumDIV[[1]])

# tidy and write out
humsumDIVTable <- resTabSimple(hum.sumDIV)
write.csv(humsumDIVTable, "data/HumSpGLMRes_sumDIV.csv", row.names = F)

mod.stk <- do.call(stack, hum.sumDIV[[3]])
writeRaster(mod.stk, file.path(mod.out.nov, "SpGLMRes_sumDIV", "hum"),format = "GTiff",
            bylayer = T, suffix = "numbers", overwrite = T)

hum.sumDIV_NoAn <- spatGLM.AnimalMod(ob.col = OB_hum_imp,
                      coV.v = c( "BB.condDIV",
                                 "BB.condDIV_2",
                                 "BB.condDIV_4",
                                 "BB.condDIV_6",
                                 "OB_ann_imp","OB_ann_imp_1", 
                                 "hdl","logPop","lnBm.div","lFrag",
                                 "OB_hum_imp","month",
                                 "x", "y", "cell"),
                      dat= dat.mod)
summary(hum.sumDIV_NoAn[[1]])
mod.stk <- do.call(stack, hum.sumDIV_NoAn[[3]])
writeRaster(mod.stk,
            file.path(mod.out.nov, "SpGLMRes_sumDIV", "humNoAn"),
            format = "GTiff",
            bylayer = T,
            suffix = "numbers",
            overwrite = T)

#### Hum qAIC ###



quack <- list(hum.full, hum.prob, hum.sum, hum.sumDIV, hum.null)
hum.qAIC <- as.data.frame(do.call(rbind, lapply(quack, qAIC)),
              row.names = c("full", "prob", "sum", "sumDIV", "null"))

write.csv(hum.qAIC, "data/humqAIC.csv")

##############################################################
#### Animal Outbreaks ####
ann.full <- spatGLM(ob.col = OB_ann_imp,
                               coV.v = c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
                                          "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                          "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                          "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                          "hdl","lnBm.div","lFrag", "month",
                                          "OB_ann_imp",  "x", "y", "cell"),
                               dat = dat)
summary(ann.full[[1]])
anFull <- resTabSimple(ann.full)
write.csv(anFull, "data/AnFullSpGLMRes_org.csv", row.names = F)

##Creating animal rasters without stocastic events
an.op <- spatGLM.AnimalMod(ob.col = OB_ann_imp, 
                           coV.v = c( "ptr_dbl_imp_BR", "mic_dbl_imp_BR", "mol_dbl_imp_BR",
                                      "ptr_dbl_imp_BR_2", "mic_dbl_imp_BR_2", "mol_dbl_imp_BR_2",
                                      "ptr_dbl_imp_BR_4", "mic_dbl_imp_BR_4", "mol_dbl_imp_BR_4",
                                      "ptr_dbl_imp_BR_6", "mic_dbl_imp_BR_6", "mol_dbl_imp_BR_6",
                                      "lnBm.div","lFrag", "month",
                                      "OB_ann_imp",  "x", "y", "cell"),
                           dat = dat)
an.stk <- do.call(stack, an.op[[3]])
writeRaster(an.stk, file.path(mod.out.nov, "SpGLMRes_ORG", "ann"),format = "GTiff",
            bylayer = T, suffix = "numbers", overwrite = T)

## Null
ann.null <- spatGLM(ob.col = OB_ann_imp,
                    coV.v = c( "lnBm.div","lFrag", "month",
                               "OB_ann_imp",  "x", "y", "cell"),
                    dat = dat)
summary(ann.null[[1]])
anNull <- resTabSimple(ann.null)
write.csv(anNull, "data/AnNullSpGLMRes_Nov.csv", row.names = F)

## alt 
ann.prob <- spatGLM(ob.col = OB_ann_imp,
                    coV.v = c( "ptr_dbl_imp", "mic_dbl_imp", "mol_dbl_imp",
                               "ptr_dbl_imp_Prob_2", "mic_dbl_imp_Prob_2", "mol_dbl_imp_Prob_2",
                               "ptr_dbl_imp_Prob_4", "mic_dbl_imp_Prob_4", "mol_dbl_imp_Prob_4",
                               "ptr_dbl_imp_Prob_6", "mic_dbl_imp_Prob_6", "mol_dbl_imp_Prob_6",
                               "hdl","lnBm.div","lFrag", "BatDiv", "month",
                               "OB_ann_imp",  "x", "y", "cell"),
                    dat = dat)
summary(ann.prob[[1]])
anAlt <- resTabSimple(ann.prob)
write.csv(anAlt, "data/AnAltSpGLMRes_prob.csv", row.names = F)

ann.sum <- spatGLM(ob.col = OB_ann_imp,
                    coV.v = c( "BB.cond",
                               "BB.cond_2",
                               "BB.cond_4",
                               "BB.cond_6",
                               "BatDiv",
                               "hdl","logPop","lnBm.div","lFrag",
                               "OB_ann_imp","month",
                               "x", "y", "cell"),
                    dat = dat)
summary(ann.sum[[1]])
anSum <- resTabSimple(ann.sum)
write.csv(anSum, "data/AnAltSpGLMRes_sum.csv", row.names = F)
## Compare the 2
fAn.qAIC <- qAIC(ann.full)
nAn.qAIC <- qAIC(ann.null)
aAn.qAIC <- qAIC(ann.prob)
sAn.qAIC <- qAIC(ann.sum)
delta.An <- nAn.qAIC - fAn.qAIC

#### Tables for publication ####
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

## human full table
hum.full.tab.pub <- humFullTable[c(5,6,4)]
hum.full.tab.pub$Significance <- c("***",
                                   "",
                                   "",
                                   "",
                                   "**",
                                   "",
                                   "",
                                   "",
                                   "",
                                   "",
                                   "*",
                                   "",
                                   "",
                                   "***",
                                   "",
                                   "",
                                   "",
                                   "***",
                                   "")
                                     
rownames(hum.full.tab.pub) <- c("Intercept",
                                "$?beta_1 ?lambda_{?afb}$",
                                "$?beta_2 ?lambda_{?mic}$",
                                "$?beta_3 ?lambda_{?mol}$",
                                "$?beta_4 ?lambda_{?afb_2}$",
                                "$?beta_5 ?lambda_{?mic_2}$",
                                "$?beta_6 ?lambda_{?mol_2}$",
                                "$?beta_7 ?lambda_{?afb_4}$",
                                "$?beta_8 ?lambda_{?mic_4}$",
                                "$?beta_9 ?lambda_{?mol_4}$",
                                "$?beta_{10} ?lambda_{?afb_6}$",
                                "$?beta_{11} ?lambda_{?mic_6}$",
                                "$?beta_{12} ?lambda_{?mol_6}$",
                                "$?beta_{13} ?mathrm{OB}_{?mathrm{an}}$",
                                "$?beta_{14} ?mathrm{OB}_{?mathrm{an}_{l-1}}$",
                                "$?beta_{15} ?mathrm{BVD}$",
                                "$?beta_{15} ?logplus(PopDen)$",
                                "$?beta_{16} ?logplus(NBM div)$",
                                "$?beta_{17} ?logplus(fragIndex)$")
  

tab <- xtable(hum.full.tab.pub)
print(tab, sanitize.text.function = function(x) {x})

## q aic table
dqaic <- as.data.frame(rbind( null.qAIC,full.qAIC, nAn.qAIC,fAn.qAIC))
dqaic$daic <- c(0,-317.23,0,-83.99)

rownames(dqaic) <- c("$RR_?mathrm{hum} Null$",
                     "$RR_?mathrm{hum} Full$",
                     "$RR_?mathrm{an} Null$",
                     "$RR_?mathrm{an} Full$")
colnames(dqaic) <- c("$k$", "$?hat{c}$", "?mathrm{qAIC}", "$?Delta ?mathrm{qAIC}$")
q.tab <- xtable(dqaic)
print(q.tab, sanitize.text.function = function(x) {x})
