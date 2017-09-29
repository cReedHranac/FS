lc <- env$LandCover
z <- Which(lc == i, cells = T)
q <- sample(z,1)
xyFromCell(lc, q)
levels(lc)

z <- list()
for( i in 1:5){
  z[[i]] <- xyFromCell(lc, sample(Which(lc ==i, cells = T), 1))
}

frame <- do.call(rbind, z)
write.csv(frame, file.path(clean.dir, "Pseudopoints.csv"), row.names = F)
