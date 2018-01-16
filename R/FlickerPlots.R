source("R/helperFunctions.R")

P.2.2 <- resRasterLoad("ptr",
                       run.id = "DBL_2",
                       dbl = T,
                       path.to.dir = mod.out.dir)
flickerPlot(P.2.2, dbl = T, path.out = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths"))

P.2.1 <- resRasterLoad("ptr",
                       run.id = "SNG_2",
                       dbl = F,
                       path.to.dir = mod.out.dir)
flickerPlot(P.2.1, dbl = T, path.out = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths_1"))




mi.2.2 <- resRasterLoad("mic",
                       run.id = "DBL_2",
                       dbl = T,
                       path.to.dir = mod.out.dir)
flickerPlot(mi.2.2, dbl = T, path.out = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths"))
mo.2.1 <- resRasterLoad("mol",
                       run.id = "SNG_2",
                       dbl = F,
                       path.to.dir = mod.out.dir)
flickerPlot(mo.2.1, dbl = T, path.out = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths_1"))


mo.2.2 <- resRasterLoad("mol",
                             run.id = "DBL_2",
                             dbl = T,
                             path.to.dir = mod.out.dir)
mi.2.1 <- resRasterLoad("mic",
                        run.id = "SNG_2",
                        dbl = F,
                        path.to.dir = mod.out.dir)
flickerPlot(mi.2.1, dbl = T, path.out = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths_1"))
flickerPlot(mi.2.2, dbl = T, path.out = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths"))










gifMaster(taxon.class = "ptr", path.to.sng = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths_1"),
          path.to.dbl = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths"), fps = .5)
gifMaster(taxon.class = "mic", path.to.sng = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths_1"),
          path.to.dbl = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths"), fps = .5)
gifMaster(taxon.class = "mol", path.to.sng = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths_1"),
          path.to.dbl = file.path("D:", "Dropbox", "FS", "FigOut", "ModeledBirths"), fps = .5)
