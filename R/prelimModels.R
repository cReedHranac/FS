#### Preliminary Modeling ####

rm(list = ls())
source("R/helperFunctions.R")
library(data.table); library(dplyr)

#### Functions ####
lagCol <- function(data.tbl, colName, lag){
  n <- nrow(data.tbl)/12
  col.v <- as.vector(data.tbl[c(seq(((12-lag)*n):nrow(data.tbl)),seq(1:((12-lag)*n-1))),colName])
  return(col.v)
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
dat <- tbl_df(fread(file.path(clean.dir, "longTable.csv")))

dat.br <- dat %>%
  mutate(ptr_BR = ptr_sng * Mega_sum,
         mic_BR = mic_sng * Micro_sum,
         mol_BR = mol_sng * Molo_sum)

dat.1 <- dat.br %>%
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
         logPop = log(popDen + 1))

#### Prelim Models ####

## AnimalModels
  #No lag
ann_mod <- glm(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR + offset(logPop), 
               data = dat.1, family = "binomial")
summary(ann_mod)

ann_mod <- glm(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR + LandCover_category + offset(logPop), 
               data = dat.1, family = "binomial")
summary(ann_mod)

  #Lag 1
ann_mod.1 <- glm(OB_ann0_ ~ ptr_BR_1 + mic_BR_1 + mol_BR_1 + offset(logPop), 
               data = dat.1, family = "binomial")
summary(ann_mod.1)

  #Lag 2
ann_mod.2 <- glm(OB_ann0_ ~ ptr_BR_2 + mic_BR_2 + mol_BR_2 + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(ann_mod.2)

  #Lag 3
ann_mod.3 <- glm(OB_ann0_ ~ ptr_BR_3 + mic_BR_3 + mol_BR_3 + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(ann_mod.3)

  #Lag 4
ann_mod.4 <- glm(OB_ann0_ ~ ptr_BR_4 + mic_BR_4 + mol_BR_4 + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(ann_mod.4)

#Lag 5
ann_mod.5 <- glm(OB_ann0_ ~ ptr_BR_5 + mic_BR_5 + mol_BR_5 + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(ann_mod.5)

#Lag 2
ann_mod.6 <- glm(OB_ann0_ ~ ptr_BR_6 + mic_BR_6 + mol_BR_6 + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(ann_mod.6)


## HumanModels
hum_mod <- glm(OB_hum0_ ~ ptr_BR + mic_BR + mol_BR + OB_ann0_/Mam_sum + offset(logPop), 
               data = dat.1, family = "binomial")
summary(hum_mod)

  #Lag 1
hum_mod.1 <- glm(OB_hum0_ ~ ptr_BR_1 + mic_BR_1 + mol_BR_1 + OB_ann0_/Mam_sum + offset(logPop), 
               data = dat.1, family = "binomial")
summary(hum_mod.1)

  #Lag 2
hum_mod.2 <- glm(OB_hum0_ ~ ptr_BR_2 + mic_BR_2 + mol_BR_2 + OB_ann0_/Mam_sum + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(hum_mod.2)

  #Lag 3
hum_mod.3 <- glm(OB_hum0_ ~ ptr_BR_3 + mic_BR_3 + mol_BR_3 + OB_ann0_/Mam_sum + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(hum_mod.3)

  #Lag 4
hum_mod.4 <- glm(OB_hum0_ ~ ptr_BR_4 + mic_BR_4 + mol_BR_4 + OB_ann0_/Mam_sum + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(hum_mod.4)

  #Lag 5
hum_mod.5 <- glm(OB_hum0_ ~ ptr_BR_5 + mic_BR_5 + mol_BR_5 + OB_ann0_/Mam_sum + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(hum_mod.5)

  #Lag 6
hum_mod.6 <- glm(OB_hum0_ ~ ptr_BR_6 + mic_BR_6 + mol_BR_6 + OB_ann0_/Mam_sum + offset(logPop), 
                 data = dat.1, family = "binomial")
summary(hum_mod.6)
