setwd("/scratch-local/gvd16/Projects_temp/2013_Priebe_Strain/WD")

source("../Code/R/PlotFunctions_for_webs.R")

base.dir <- "../Data/Elaborated/interactions"


for(folder in base.dir %>% list.dirs){

folder %>%
  {paste0("../Plots/",.,".pdf")} %>%
  pdf(width=9,height=3,useDingbats=FALSE)
  
c(0, 0, 0, 0) %>% {par(mar= . )}

  base.dir %>% 
    paste0("/", folder, "/interactions.csv") %>% #substitute with 'folder' for loooping
    read.csv(head=FALSE) %>%
    .[,c(2,1)] %>%
    graph.data.frame(directed = TRUE) %>%
    simplify(remove.loops = FALSE) %>%
    Plot.foodweb(dodge=0.1, size.scale=2)


dev.off()

}
