#### fitting for spastat

source("R/helperFunctions.R")
library(data.table); library(dplyr)

## Dave run this!
clean.dir <- "."


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
# te <- tess(tiles=windows)


library(dplyr)
dat.sel <- dat.2 #%>%
  #filter( OB_hum0_ == 1)


coordinates(dat.sel) <- ~ x + y
proj4string(dat.sel) <- CRS(wgs)

ob.pts <- dat.2  %>%
  filter( OB_hum0_ == 1)
coordinates(ob.pts) <- ~ x + y
proj4string(ob.pts) <- CRS(wgs)

## Need to filter with unique points
ob.pts <- remove.duplicates(ob.pts)

ob <- ppp(ob.pts@coords[,"x"],ob.pts@coords[,"y"],
          window = windows[[1]])
# plot(ob)
# summary(ob)

### STUFF FOR DETERMINING WEIGHTS etc
im_win <- as.im(rf)
quad_scheme <- pixelquad(ob, im_win)
p_for_weights <- ppm(quad_scheme)

# code under here is basically doing: ppm(quad_scheme) except without running ppm (i.e.
# it's doing a GLM with the appropriate weights using dummy points etc)

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

# Adding covariates: Just left join on x,y

dat.list <- dat.2[, c( "ptr_BR", "mic_BR", "mol_BR",
                        "ptr_BR_1", "mic_BR_1", "mol_BR_1",
                        "ptr_BR_3", "mic_BR_3", "mol_BR_3",
                        "ptr_BR_5", "mic_BR_5", "mol_BR_5",
                        "logPop", "OB_ann0_", "NB_lDiv","month",
                       "OB_hum0_",  "x", "y")]
glm.full <- left_join(glmdata, dat.list, by = c("x", "y"))

# Filtering for month, and adjusting responce weight
glm.m <- glm.full %>%
  filter(month == 12) %>%
  mutate(.mpl.Y = .mpl.Y*OB_hum0_) #using binary responce vector to clear data from some points

#for month sub sets
glmdata <- glm.m
#for full data
glmdata <- glm.full

.mpl.W <- glmdata$.mpl.W
.mpl.SUBSET <- NULL

# This is where we'd put covariaty stuff in somehow
# rhs <- '~x + y'
# fmla <- paste(".mpl.Y ", rhs)
# fmla <- as.formula(fmla)
int.form <- as.formula(paste(".mpl.Y ","~", paste(names(dat.list[1:(length(names(dat.list))-4)]), collapse = "+")))

fmla <- int.form

# and running the model fit
FIT  <- glm(fmla, family=quasi(link="log", variance="mu"), weights=.mpl.W, data=glmdata)

