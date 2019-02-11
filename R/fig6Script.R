source("R/helperFunctions.R")
library(ggplot2); library(dplyr); library(data.table); library(gridExtra); library(gtools)
library(rgdal); library(raster); library(ggridges); library(RcppRoll)
library(latex2exp); library(patchwork)

#### functions ####
better.names <- function(x, model.name){
  ### function for impoving names accociated with items retrieved from SpatHandler
  base <- sapply(strsplit(model.name, "_"), tail, 1)
  i <- 1:12
  j <- c(i[12],i[1:11])
  names(x) <- paste0(base, "_", j, "_", i)
  return(x)
}

spatHandler <- function(model.name, mod.dir){
  ## function for loading rasters produced from spatGLM and producing an averaged product based on the
  ## model string argument
  base <- sapply(strsplit(model.name, "_"), tail, 1)
  
  f.list <- mixedsort(list.files(file.path(mod.dir,model.name),
                                 pattern = "noAn_2_", 
                                 full.names = T))
  
  if(!any(file.exists(f.list))){
    stop("string not found try again \n *cough* dumbass *cough*")
  }
  ## order and read
  stk <- better.names(stack(f.list), model.name = model.name)
  m.stk <- mean(stk)
  out.l <- list(stk,m.stk)
  return(out.l)
}

model_names <- tribble(~Mod.Name, ~Plot.Name,
        'null',  'NULL',
        'nDiv', 'D_{tot}',
        'cBDiv', 'D+B_{tot}',
        'cNDiv', 'B_{tot}',
        'cPDiv', 'D \\times B_{tot}',
        'nSBD', 'D_{tax}',
        'Prb', 'D + B_{tax}',
        'ORG', 'D \\times B_{tax}'
        )
model_names$Plot.Name = factor(model_names$Plot.Name, levels=model_names$Plot.Name)

altModBoxes.Month <- function(x, df = ob.masterframe,
                              write.out = F, out.fig = NULL){
  ### Function for looking at outbreak prediction preformace across alternative models
  ## x <- name of colum from df to use
  ## df <- df to use (default to ob.masterframe)
  ## write.out <- logical, should write
  ## out.fig <- where the figure should be written to
  
  x.enquo <- enquo(x)
  x.quo <- quo(x)
  
  ob.place.time <- df %>%
    filter(Outbreak == "Beni" & window == 7)
  ob.pt <- df %>% 
    filter(Outbreak == "Bikoro" & window == 4) %>%
    bind_rows(ob.place.time)
  
  p.vol <- ggplot(ob.pt,
                  aes(x = Outbreak,
                      y = !!x.enquo,
                      color = Mod.Name)) + 
    geom_boxplot(size = 1, width = .5) + 
    # geom_jitter(shape = 16, 
    #             position = position_jitter(.1),
    #             color = "black")+
    #Need to figure out why this doesnt work
    # scale_y_continuous(limits = c(floor(!!x.enquo),
    #                               ceiling(!!x.quo))) + 
    labs(x = "Outbreak", 
         y = quo_name(x.enquo))+
    theme_bw()+
    theme(#legend.position = "none",
      axis.title.x = element_blank())
  
  if(write.out == T){
    a <- file.path(out.fig, paste0("altModBoxesMonth_",quo_name(x.enquo),".pdf"))
    ggsave(a,
           plot = p.vol,
           device = cairo_pdf,
           height = 7,
           width = 7.5,
           units = "in",
           dpi = 300)
  }
  return(p.vol)
  
}
altModBoxes.Total <- function(x, df = ob.masterframe,
                              write.out = F, out.fig = NULL){
  ### Function for looking at outbreak prediction preformace across alternative models
  ## x <- name of colum from df to use
  ## df <- df to use (default to ob.masterframe)
  ## write.out <- logical, should write
  x.enquo <- enquo(x)
  x.quo <- quo(x)
  
  p.vol <- ggplot(df,
                  aes(x = Outbreak,
                      y = !!x.enquo,
                      color = Mod.Name)) + 
    geom_boxplot(size = 1, width = .5) + 
    # geom_jitter(shape = 16, 
    #             position = position_jitter(.1),
    #             color = "black")+
    #Need to figure out why this doesnt work
    # scale_y_continuous(limits = c(floor(!!x.enquo),
    #                               ceiling(!!x.quo))) + 
    labs(x = "Outbreak", 
         y = quo_name(x.enquo))+
    theme_bw()+
    theme(#legend.position = "none",
      axis.title.x = element_blank())
  
  if(write.out == T){
    a <- file.path(out.fig, paste0("altModBoxesTotal_",quo_name(x.enquo),".pdf"))
    ggsave(a,
           plot = p.vol,
           device = cairo_pdf,
           height = 7,
           width = 7.5,
           units = "in",
           dpi = 300)
  }
  return(p.vol)
  
}

#### ####
human.model.names <- c("h_cBDiv", "h_cNDiv", "h_cPDiv", "h_nDiv",
                       "h_nSBD",  "h_null",  "h_ORG",   "h_Prb"  )
## read in the dataframes in and add Model name column
## get names 
dfs <- lapply(list.files(dOut.1, pattern = "ob", full.names = T),read.csv)


for(i in 1:length(human.model.names)){
  dfs[[i]]$Mod.Name <- sapply(strsplit(human.model.names[[i]], "_"), tail, 1)
}

ob.masterframe <- do.call(rbind, dfs)
names(ob.masterframe)

t1 <- altModBoxes.Month(pct.rank,
                        write.out = F)
t2 <- altModBoxes.Month(rel.Risk,
                        write.out = F)

t3 <- altModBoxes.Total(pct.rank,
                        write.out = F)
t4 <- altModBoxes.Total(rel.Risk,
                        write.out = F)

# filter resuls down to month data
ob.month <- ob.masterframe %>% filter((Outbreak == "Beni" & window == 7) |
                          (Outbreak == "Bikoro" & window == 4)) %>% left_join(model_names)

g1 = ggplot(ob.month, aes(x = Outbreak, y = pct.rank, fill = Plot.Name)) +
  geom_boxplot(width = 0.5, position = position_dodge(width=0.7)) +
  labs(x = "Outbreak", 
       y = "Percent rank") +
  theme_bw()+
  theme(#legend.position = "none",
    axis.title.x = element_blank()) +
  guides(fill = 'none')

g2 = ggplot(ob.month, aes(x = Outbreak, y = rel.Risk, fill = Plot.Name)) +
  geom_boxplot(width = 0.5, position = position_dodge(width=0.7)) +
  labs(x = "Outbreak", 
       y = "Relative risk") +
  theme_bw()+
  theme(#legend.position = "none",
    axis.title.x = element_blank(),
      legend.text.align = 0) +
  scale_fill_discrete(name="Model", labels = unname(TeX(unique(paste0("$",model_names$Plot.Name)))))

fig6 <- g1 + g2

ggsave(filename = file.path(fig.pub, "Fig6.pdf"), 
       plot = fig6,
       device = cairo_pdf,
       width = 8,
       height = 4,
       units = "in",
       dpi = 300)