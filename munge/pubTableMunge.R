#### Publication Tables #### 
#### Tables for publication ####
source("R/helperFunctions.R")
library(tidyverse);
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

main <- read.csv(file.path(dOut.1, "top2Table.csv"))

main.sub <- main[,c(5,6,4, 10,11,9)]

rownames(main.sub) <- c("Intercept",
                     "$?beta_1 ?logplus(P_{?afb})$",
                     "$?beta_2 ?logplus(P_{?mic})$",
                     "$?beta_3 ?logplus(P_{?mol})$",
                     "$?beta_4 ?logplus(P_{?afb_2})$",
                     "$?beta_5 ?logplus(P_{?mic_2})$",
                     "$?beta_6 ?logplus(P_{?mol_2})$",
                     "$?beta_7 ?logplus(P_{?afb_4})$",
                     "$?beta_8 ?logplus(P_{?mic_4})$",
                     "$?beta_9 ?logplus(P_{?mol_4})$",
                     "$?beta_{10} ?logplus(P_{?afb_6})$",
                     "$?beta_{11} ?logplus(P_{?mic_6})$",
                     "$?beta_{12} ?logplus(P_{?mol_6})$",
                     "$?beta_{13} ?logplus(Div_{?afb})$",
                     "$?beta_{14} ?logplus(Div_{?mic})$",
                     "$?beta_{15} ?logplus(Div_{?mol})$",
                     "$?beta_{16} ?logplus(Div_{?mathrm{NBM}})$",
                     "$?beta_{17} ?logplus(PopDen)$",
                     "$?beta_{18} ?logplus(fragIndex)$",
                     "$?beta_{19} ?mathrm{BVD}$",
                     "$?beta_{20} ?mathrm{OB}_{?mathrm{an}}$",
                     "$?beta_{21} ?mathrm{OB}_{?mathrm{an}_{l-1}}$")

tab <- xtable(main.sub)
print(tab, sanitize.text.function = function(x) {x})



