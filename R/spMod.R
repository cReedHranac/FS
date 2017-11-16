#### fitting for spastat

rm(list = ls())
source("R/helperFunctions.R")
library(data.table); library(dplyr)

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



#### Data ####
dat <- tbl_df(fread(file.path(clean.dir, "longTable.csv")))

dat.br <- dat %>%
  mutate(ptr_BR = log(ptr_dbl * Mega_sum +1),
         mic_BR = log(mic_dbl * Micro_sum +1),
         mol_BR = log(mol_dbl * Molo_sum +1))

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
         mol_BR_6 = wrap(mol_BR, n=6, order_by = month)) %>%
  ungroup %>%
  mutate(logPop = log(popDen + 1),
         NB_lDiv = log((Mam_sum - (Mega_sum + Molo_sum + Micro_sum))+1))


## Create Owin object ##
library(raster) ; library(maptools) ; library(data.table)
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rf <- raster(file.path(data.source, "CropMask.tif"))
rf.poly <- rasterToPolygons(rf, fun = function(x){x==1}, dissolve = T)

library(spatstat)
regions <- slot(rf.poly, "polygons")
regions <- lapply(regions, function(x) { SpatialPolygons(list(x)) })
windows <- lapply(regions, as.owin)
te <- tess(tiles=windows)

## lets try 
library(dplyr)
dat.sel <- dat.2 #%>%
  #filter( OB_hum0_ == 1)


coordinates(dat.sel) <- ~ x + y
proj4string(dat.sel) <- CRS(wgs)

ob.pts <- dat.2  %>%
  filter( OB_hum0_ == 1)
coordinates(ob.pts) <- ~ x + y
proj4string(ob.pts) <- CRS(wgs)


ob <- ppp(ob.pts@coords[,"x"],ob.pts@coords[,"y"],
          window = windows[[1]])
plot(ob)
summary(ob)

dat.list <- dat.sel[, c("ptr_BR", "mic_BR", "mol_BR",
                      "ptr_BR_1", "mic_BR_1", "mol_BR_1",
                      "ptr_BR_3", "mic_BR_3", "mol_BR_3",
                      "ptr_BR_5", "mic_BR_5", "mol_BR_5",
                      "logPop", "OB_ann0_", "NB_lDiv")]
hf <- hyperframe(window = windows)
hf <- cbind.hyperframe(hf, dat.list@data)

int.form <- as.formula(paste("OB_hum0_", "~", paste(names(dat.list[2:length(names(dat.list))]), collapse = "+")))

ft.int <- ppm(ob, trend = as.formula(int.form), covariates = hf)
