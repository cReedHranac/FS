#### MIC 9 Attemtpts####
source("R/helperFunctions.R")
source("R/bioMod_fun.R")

occLoad("mic")

# set.seed(5) That makes it fail first model run

## if I just let it run normally then it 1 time miraciously ran all the way though
bioMod(mod.id = "mic9", occ.db = mic.occ, nrep = 10, run.id = "SNG_2", dbl = F)
