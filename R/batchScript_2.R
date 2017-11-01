source("R/helperFunctions.R")
classMod("ptr", run.id = "DBL_2", nrep = 10, dbl = T)
classMod("mic", run.id = "DBL_2", nrep = 10, dbl = T)
classMod("mol", run.id = "DBL_2", nrep = 10, dbl = T)

classMod("ptr", run.id = "SNG_2", nrep = 10, dbl = F)
classMod("mic", run.id = "SNG_2", nrep = 10, dbl = F)
classMod("mol", run.id = "SNG_2", nrep = 10, dbl = F)

source("R/dblMonthFuns.R")
source("R/bioMod_fun.R")
source("R/helperFunctions.R")


mic.occ <- occLoad2("mic")
mic.l <- listGen(mic.occ)
lapply(mic.l, bioMod, occ.db=mic.occ, nrep=10, run.id = "DBL_2", dbl=T)
rm(list=ls())

source("R/dblMonthFuns.R")
source("R/bioMod_fun.R")
source("R/helperFunctions.R")
occLoad("ptr")
ptr.l <- listGen(ptr.occ)
lapply(ptr.l, bioMod, occ.db=ptr.occ, nrep=10, run.id = "SNG_2", dbl=F)
rm(list=ls())

source("R/dblMonthFuns.R")
source("R/bioMod_fun.R")
source("R/helperFunctions.R")
occLoad("mic")
mic.l <- listGen(mic.occ)
lapply(mic.l, bioMod, occ.db=mic.occ, nrep=10, run.id = "SNG_2", dbl=F)
bioMod("mic10",occ.db=mic.occ, nrep=10, run.id = "SNG_2", dbl=F)
bioMod("mic11",occ.db=mic.occ, nrep=10, run.id = "SNG_2", dbl=F)
bioMod("mic12",occ.db=mic.occ, nrep=10, run.id = "SNG_2", dbl=F)
  ## causes R session to abort and crash
# bioMod("mic9",occ.db=mic.occ, nrep=10, run.id = "SNG_2", dbl=F)
rm(list=ls())


source("R/bioMod_fun.R")
source("R/helperFunctions.R")
occLoad("mol")
mol.l <- listGen(mol.occ)
lapply(mol.l, bioMod, occ.db=mol.occ, nrep=10, run.id = "SNG_2", dbl=F)
rm(list=ls())
