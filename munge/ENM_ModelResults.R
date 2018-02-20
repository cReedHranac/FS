#### Extracting and forming the results from biomod modeling ####
source("R/helperFunctions.R")
library(dplyr)
mod.path <- file.path(mod.out.dir, "DBL_2")
m.list <- list.dirs(mod.path, recursive = F, full.names = F)
m.path <- list.dirs(mod.path, recursive = F)

res.list <- list()
for(i in 1:length(m.list)){
  dat <- read.csv(file.path(m.path[[i]], paste0(m.list[[i]], "_EM_Eval.csv")))
  sel <- dat[12,]
  sel$m <- m.list[[i]]
  res.list[[i]] <- sel
}
res <- do.call(rbind, res.list)

z <- data.frame(Model = res$m,
                ROC = res$Testing.data,
                Sensitivity = res$Sensitivity,
                Specificity = res$Specificity)
library(gtools)
ord <- mixedorder(m.list)
z.1 <- z[ord,]
write.csv(z.1, file.path(clean.dir, "EMN_ModelResults.csv"), row.names = F)
