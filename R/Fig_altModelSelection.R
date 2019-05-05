source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)
library(latex2exp); library(patchwork)

theme_set(theme_bw(base_family = "serif"))


data_outdf   <- mod.out.nov

# model name formatting
model_names <- tribble(~Mod.Name, ~Plot.Name,
        'null',  'NULL',
        'nDiv', 'D_{tot}',
        'nSBD', 'D_{tax}',
        'cBDiv', 'D + B_{tot}',
        'cNDiv', 'B_{tot}',
        'cPDiv', 'D \\times B_{tot}',
        'mcBDiv', 'D + B_{tot}*',
        'mcNDiv', 'B_{tot}*',
        'mcPDiv', 'D \\times B_{tot}*',
        'ORG', 'D \\times B_{tax}',
        'mORG', 'D \\times B_{tax}*',
        'Prb', 'D + B_{tax}',
        'mod', 'D + B_{tax}*')
model_names$Plot.Name = factor(model_names$Plot.Name, levels=model_names$Plot.Name)

#### ####
human.model.names <-  c("h_null", "h_nDiv","h_nSBD",
                        "h_cBDiv", "h_cNDiv", "h_cPDiv",
                        "h_mcBDiv", "h_mcNDiv", "h_mcPDiv",
                        "h_ORG","h_mORG",  "h_Prb", "h_Mod"  )
## read in the dataframes in and add Model name column
## get names 
dfs <- lapply(paste0(data_outdf, 
                     "/obDF", human.model.names, ".csv"),
              read.csv)


for(i in 1:length(human.model.names)){
   dfs[[i]]$Mod.Name <- sapply(strsplit(human.model.names[[i]], "_"), tail, 1)
}

ob.masterframe <- do.call(rbind, dfs)
names(ob.masterframe)

# filter resuls down to month data
ob.month <- ob.masterframe %>% filter((Outbreak == "Beni" & window == 7) |
                          (Outbreak == "Bikoro" & window == 4)) %>% left_join(model_names)

g1 = ggplot(ob.month, aes(x = Outbreak, y = pct.rank, fill = Plot.Name)) +
  geom_boxplot(width = 0.5, position = position_dodge(width=0.7)) +
  labs(x = "Outbreak", 
       y = "Percent rank") +
  theme(#legend.position = "none",
    axis.title.x = element_blank()) +
  guides(fill = 'none')

g2 = ggplot(ob.month, aes(x = Outbreak, y = rel.Risk, fill = Plot.Name)) +
  geom_boxplot(width = 0.5, position = position_dodge(width=0.7)) +
  labs(x = "Outbreak", 
       y = "Relative risk") +
  scale_y_log10() +
  theme(#legend.position = "none",
    axis.title.x = element_blank(),
      legend.text.align = 0) +
  scale_fill_discrete(name="Model", labels = unname(TeX(unique(paste0("$",model_names$Plot.Name)))))

(fig6 <- g1 + g2)

ggsave(filename = file.path(fig.pub, "Fig_altModelSeltection.pdf"), 
       plot = fig6,
       device = cairo_pdf,
       width = 8,
       height = 4,
       units = "in",
       dpi = 300)
