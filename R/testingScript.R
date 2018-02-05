occLoad("mic")
apply(mic.occ[3:ncol(mic.occ)],2, sum)

bz <- long.table.br %>% dplyr::group_by(cell) %>%
  dplyr::mutate(mol_dbl_imp_BR_1 := wrap(mol_dbl_imp_BR,  n= l, order_by = month)) %>%
  ungroup %>%
dplyr::select(mol_dbl_imp_BR_1)
bz
rm (list = ls())

library(dplyr)

test <- data.frame(t1=rep(1:3, each=3), t2=rep(letters[1:3], 3))
var_name1 <- "t1"
var_name2 <- "t2"
var_name3 <- "comb"

foo <- list(lazyeval::interp(~paste(x, y), x=as.name(var_name1), y=as.name(var_name2)))
names(foo) <- var_name3


# one liner FTW!
test %>% mutate_(.dots = setNames(list(lazyeval::interp(~paste(x, y), x=as.name(var_name1), y=as.name(var_name2))), var_name3))

z <- spatGLM(ob.col = "OB_hum_imp",
             coV.v = c( "ptr_dbl_raw_BR", "mic_dbl_raw_BR", "mol_dbl_raw_BR",
                                            "ptr_dbl_raw_BR_2", "mic_dbl_raw_BR_2", "mol_dbl_raw_BR_2",
                                            "ptr_dbl_raw_BR_4", "mic_dbl_raw_BR_4", "mol_dbl_raw_BR_4",
                                            "ptr_dbl_raw_BR_6", "mic_dbl_raw_BR_6", "mol_dbl_raw_BR_6",
                                            "logPop", "OB_ann_imp", "NB_lDiv",
                                            "fragIndex", "OB_ann_imp_1","month",
                                            "OB_hum_imp",  "x", "y", "cell"),
             dat = dat, 
             rGrid = rf)

ptr.occ<- occLoad2("ptr")
mic.occ <- occLoad2("mic")
mol.occ <- occLoad2("mol")

z.1 <-apply(ptr.occ[3:dim(ptr.occ)[[2]]], 2, sum)
z.2 <-apply(mic.occ[3:dim(mic.occ)[[2]]], 2, sum)
z.3 <-apply(mol.occ[3:dim(mol.occ)[[2]]], 2, sum)
q <- rbind(z.1, z.2, z.3)
rownames(q) <- c("ptr", "mic", "mol")
colnames(q) <- 1:12
write.csv(z, file = "D://Dropbox/reed/PrelimModels/birthSummaryTable.csv")
