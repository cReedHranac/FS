#### Pseudo Point Generator ####
source("R/helperFunctions.R")

env <- staticStkLoad()

lc <- env$LandCover

## For psedudo occurences
z <- list()
for( i in 1:5){
  z[[i]] <- xyFromCell(lc, sample(Which(lc ==i, cells = T), 1))
}
## For pseudo absences
q <- list()
for( i in 1:5){
  q[[i]] <- xyFromCell(lc, sample(Which(lc ==i, cells = T), 1))
}
## 5 more?
w <- list()
for( i in 1:5){
  w[[i]] <- xyFromCell(lc, sample(Which(lc ==i, cells = T), 1))
}
frame <- do.call(rbind, c(z,q,w))
write.csv(frame, file.path(clean.dir, "Pseudopoints.csv"), row.names = F)
