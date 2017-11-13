#### Preliminary Modeling ####

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

dat.1 <- dat.br %>% group_by(cell) %>%
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
dat.1 <- as.data.frame(dat.1)
#### Prelim Models ####

## AnimalModels
  #No lag 0 back  
# ann_mod <- glm(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR + logPop, 
#                data = dat.1, family = binomial(link = "log"))
# summary(ann_mod)
# 
# an_mod_full <- glm(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR +
#                      ptr_BR_1 + mic_BR_1 + mol_BR_1 +
#                      ptr_BR_2 + mic_BR_2 + mol_BR_2 +
#                      ptr_BR_3 + mic_BR_3 + mol_BR_3 +
#                      ptr_BR_4 + mic_BR_4 + mol_BR_4 +
#                      ptr_BR_5 + mic_BR_5 + mol_BR_5 +
#                      ptr_BR_6 + mic_BR_6 + mol_BR_6 +
#                      logPop,data= dat.1, family = binomial(link = "log"))
# 
# summary(an_mod_full)

  ####GLMNET
library(glmnet)
x <- model.matrix(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR +
                    ptr_BR_1 + mic_BR_1 + mol_BR_1 +
                    ptr_BR_2 + mic_BR_2 + mol_BR_2 +
                    ptr_BR_3 + mic_BR_3 + mol_BR_3 +
                    ptr_BR_4 + mic_BR_4 + mol_BR_4 +
                    ptr_BR_5 + mic_BR_5 + mol_BR_5 +
                    ptr_BR_6 + mic_BR_6 + mol_BR_6 +
                    logPop,data= dat.1)
y <- as.matrix(dat.1[as.numeric(rownames(x)),"OB_ann0_"])


cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 5) # up to 50
plot(cvfit)
cvfit$lambda.1se
## ^^ WONT RUN NOW


x <- model.matrix(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR +
                    logPop,data= dat.1)
y <- as.matrix(dat.1[as.numeric(rownames(x)),"OB_ann0_"])

cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 5)
plot(cvfit)

## HumanModels
# hum_mod <- glm(OB_hum0_ ~ ptr_BR + mic_BR + mol_BR + OB_ann0_/Mam_sum +logPop, 
#                data = dat.1, family = binomial(link = "log"))
# summary(hum_mod)

x <- model.matrix(OB_hum0_ ~ ptr_BR + mic_BR + mol_BR +
                  OB_ann0_/Mam_sum +logPop, 
                  data = dat.1)
y <- as.matrix(dat.1[as.numeric(rownames(x)),"OB_hum0_"])
cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 5)
plot(cvfit)


x <- model.matrix(OB_hum0_ ~ ptr_BR + mic_BR + mol_BR +
                    ptr_BR_1 + mic_BR_1 + mol_BR_1 +
                    ptr_BR_2 + mic_BR_2 + mol_BR_2 +
                    ptr_BR_3 + mic_BR_3 + mol_BR_3 +
                    ptr_BR_4 + mic_BR_4 + mol_BR_4 +
                    OB_ann0_/Mam_sum, + logPop, 
                  data = dat.1)
y <- as.matrix(dat.1[as.numeric(rownames(x)),"OB_hum0_"])
hum.fit <- glmnet(x, y, family = "binomial")
plot(hum.fit)



