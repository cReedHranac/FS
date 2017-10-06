# Run code for *_1 and the figure production

source("R/helperFunctions.R")
# classMod("ptr", run.id = "DBL_1", nrep = 10, dbl = T)
# classMod("mic", run.id = "DBL_1", nrep = 10, dbl = T)
classMod("mol", run.id = "DBL_1", nrep = 10, dbl = T)

classMod("ptr", run.id = "SNG_1", nrep = 10, dbl = F)
classMod("mic", run.id = "SNG_1", nrep = 10, dbl = F)

classMod(taxon.class = "mol", run.id = "SNG_1", nrep = 10, dbl = F)


DBL_1 <- file.path("D:", "Dropbox", "FS", "Figout", "DBL_1")
SNG_1 <- file.path("D:", "Dropbox", "FS", "Figout", "SNG_1")

ptr2.stk <- resRasterLoad("ptr", "DBL_1", T, "ModOut")
flickerPlot(ptr2.stk, T, births = T,
                          path.out = file.path(DBL_1, "PTR_Birth"))
# 
# mic2.stk <- resRasterLoad("mic", "DBL_1", T, "ModOut")
# flickerPlot(mic2.stk, T, births = T, 
#             path.out = file.path(DBL_1, "MIC_Birth"))
# 
# mol2.stk <- resRasterLoad(taxon.class = "mol", run.id = "DBL_1",dbl =  T, path.to.dir = "ModOut")
# ## modify txcl to get the correct list vector (removing 12.1 for now since the model fails)
# txc.l <- txc.l[-12]
# mol2.stk <- r.stk
# 
# flickerPlot(mol2.stk, T, births = T,
#             path.out = file.path(DBL_1, "MOL_Birth"))
# 

# 
ptr1.stk <- resRasterLoad("ptr", "SNG_1", F, "ModOut")
flickerPlot(ptr1.stk, F, births = T,
            path.out = file.path(SNG_1, "PTR_Birth"))
# 
# mic1.stk <- resRasterLoad("mic", "SNG_1", F, "ModOut")
# flickerPlot(mic1.stk, F, births = T, 
#             path.out = file.path(SNG_1, "MIC_Birth"))
# 
# mol1.stk <- resRasterLoad("mol", "SNG_1", F, "ModOut")
# flickerPlot(mol1.stk, F, births = T, 
#             path.out = file.path(SNG_1, "MOL_Birth"))
# 



### Fails at mol12.mol1, build prelim models 
bioMod(mod.id = z[[12]],occ.db =  occLoad2("mol"), nrep = 10, run.id = "DBL_1", dbl = T)
# Warning in .Biomod.Models.check(Model, Data, Options, calibLines, Yweights,  :
#                                   mol12.mol1 GLM was switch off because of no both
#Error in checkForRemoteErrors(val) : 
 # 3 nodes produced errors; first error: The dataset size is too small or subsampling rate is too large:
#nTrain*bag.fraction <= n.minobsinnode
#Observed or fited data contains a unique value.. Be carefull with this models predictions

DBL_1 <- file.path("D:", "Dropbox", "FS", "Figout", "DBL_1")
SNG_1 <- file.path("D:", "Dropbox", "FS", "Figout", "SNG_1")
ot <- file.path("D:", "Dropbox", "FS", "Figout", "gif_1")
## PTR
gifMaster(taxon.class = "ptr", path.to.dbl = file.path(DBL_1, "PTR_Birth"),fps = 1, master = F,
          path.out = ot)
gifMaster(taxon.class = "ptr", path.to.sng = file.path(SNG_1, "PTR_Birth"),fps = 1, master = F,
               path.out = ot)
gifMaster(taxon.class = "ptr", path.to.sng = file.path(SNG_1, "PTR_Birth"), path.to.dbl = file.path(DBL_1, "PTR_Birth"),
          fps = 1, master = T,
          path.out = ot)
## MIC
gifMaster(taxon.class = "mic", path.to.dbl = file.path(DBL_1, "MIC_Birth"),fps = 1, master = F,
          path.out = ot)
gifMaster(taxon.class = "mic", path.to.sng = file.path(SNG_1, "MIC_Birth"),fps = 1, master = F,
          path.out = ot)
gifMaster(taxon.class = "mic", path.to.sng = file.path(SNG_1, "MIC_Birth"), path.to.dbl = file.path(DBL_1, "MIC_Birth"),
          fps = 1, master = T,
          path.out = ot)
## MOl
gifMaster(taxon.class = "mol", path.to.dbl = file.path(DBL_1, "MOL_Birth"),fps = 1, master = F,
          path.out = ot)
gifMaster(taxon.class = "mol", path.to.sng = file.path(SNG_1, "MOL_Birth"),fps = 1, master = F,
          path.out = ot)
gifMaster(taxon.class = "mol", path.to.sng = file.path(SNG_1, "MOL_Birth"), path.to.dbl = file.path(DBL_1, "MOL_Birth"),
          fps = 1, master = T,
          path.out = ot)