#### Preliminary Modeling ####

rm(list = ls())
source("R/helperFunctions.R")
library(data.table); library(dplyr)

#### Data ####
dat <- tbl_df(fread(file.path(clean.dir, "longTable.csv")))

dat.br <- dat %>%
  mutate(ptr_BR = ptr_sng * Mega_sum,
         mic_BR = mic_sng * Micro_sum,
         mol_BR = mol_sng * Molo_sum)



nrpm <- nrow(dat)/12 # number rows per month
dat.lag4 <- dat.br
# maybe check the sequence hack... mightbe off by one...
dat.lag4$ptr_BR4 <- c(dat.lag4$ptr_BR[seq((8*nrpm):nrow(dat.lag4))], dat.lag4$ptr_BR[seq(1:(8*nrpm-1))])
dat.lag4$mic_BR4 <- c(dat.lag4$mic_BR[seq((8*nrpm):nrow(dat.lag4))], dat.lag4$mic_BR[seq(1:(8*nrpm-1))])
dat.lag4$mol_BR4 <- c(dat.lag4$mol_BR[seq((8*nrpm):nrow(dat.lag4))], dat.lag4$mol_BR[seq(1:(8*nrpm-1))])

dat.lag3 <- dat.br
# maybe check the sequence hack... mightbe off by one...
dat.lag3$ptr_BR3 <- c(dat.lag3$ptr_BR[seq((9*nrpm):nrow(dat.lag3))], dat.lag3$ptr_BR[seq(1:(9*nrpm-1))])
dat.lag3$mic_BR3 <- c(dat.lag3$mic_BR[seq((9*nrpm):nrow(dat.lag3))], dat.lag3$mic_BR[seq(1:(9*nrpm-1))])
dat.lag3$mol_BR3 <- c(dat.lag3$mol_BR[seq((9*nrpm):nrow(dat.lag3))], dat.lag3$mol_BR[seq(1:(9*nrpm-1))])

#### Prelim Models ####

## AnimalModels
ann_mod <- glm(OB_ann0_ ~ ptr_BR + mic_BR + mol_BR + OB_ann0_/log(Mam_sum) + offset(log(popDen)), 
               data = dat.br, family = "binomial")
summary(ann_mod)

##PTR only
ann_modPTR <- glm(OB_ann0_ ~ ptr_BR + OB_ann_bi/log(Mam_sum) + offset(log(popDen)), 
               data = dat.br, family = "binomial")
summary(ann_modPTR)


## HumanModels
hum_mod <- glm(OB_hum0_ ~ ptr_BR + mic_BR + mol_BR + OB_ann0_/log(Mam_sum) + offset(log(popDen)), 
               data = dat.br, family = "binomial")
summary(hum_mod)

  ## 4 month lag
hum_mod4 <- glm(OB_hum ~ ptr_BR4 + mic_BR4 + mol_BR4 + OB_ann_bi/log(Mam_sum) + offset(log(popDen)), 
               data = dat.lag4, family = "binomial")
summary(hum_mod4)

  ## 3 month lag
hum_mod3 <- glm(OB_hum ~ ptr_BR3 + mic_BR3 + mol_BR3 + OB_ann_bi/log(Mam_sum) + offset(log(popDen)), 
                data = dat.lag3, family = "binomial")
summary(hum_mod3)





