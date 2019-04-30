### batchScript3 
# Third times the charm?
source("R/helperFunctions.R")
classMod("ptr", run.id = "DBL_3", nrep = 10, dbl = T)
# classMod("mic", run.id = "DBL_3", nrep = 10, dbl = T)
# classMod("mol", run.id = "DBL_3", nrep = 10, dbl = T)
# 
# ##Need to find the right seed for this one
# # literally just keep attempting to run it untill it works :)
# source("R/dblMonthFuns.R"); source("R/bioMod_fun.R")
# z <- occLoad2("mic")
# bioMod("mic8.mic9", z, nrep = 10, run.id = "DBL_3")
# 
# source("R/helperFunctions.R")
# classMod("ptr", run.id = "SNG_3", nrep = 10, dbl = F)
# classMod("mic", run.id = "SNG_3", nrep = 10, dbl = F)
# classMod("mol", run.id = "SNG_3", nrep = 10, dbl = F)
# 
# ##Need to find the right seed for this one
# # literally just keep attempting to run it untill it works :)
# source("R/helperFunctions.R")
# source("R/bioMod_fun.R")
# occLoad("mic")
# bioMod("mic9", mic.occ, nrep = 10, run.id = "SNG_3", dbl = F)
# 
