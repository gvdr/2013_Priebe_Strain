"/scratch-local/gvd16/Projects_temp/2013_Priebe_Strain/WD" %>% setwd

base.dir <- "../Data/Elaborated/interactions"

"../Code/R/PlotFunctions_for_webs.R" %>% source

for(folder in base.dir %>% list.dirs){
  
  folder %>%
  {paste0("../Plots/",.,"distr.pdf")} %>%
  pdf(width=9,height=3,useDingbats=FALSE)

  base.dir %>% 
    paste0("/", folder, "/interactions.csv") %>% #substitute with 'folder' for loooping
    read.csv(head=FALSE) %>%
    .[,c(2,1)] %>%
    graph.data.frame(directed = TRUE) %>%
    simplify(remove.loops = FALSE) %>%
    get.adjacency(type="both",sparse = FALSE) %>%
    as.matrix %>%
    {
      data.frame(
      Graph = folder,
      Strains.Total = compute_strain(.,3, mode = "total"),
      Strains.Inward = compute_strain(.,3, mode = "inward"),
      Strains.Outward = compute_strain(.,3, mode = "outward")
      )
    } %>%
    melt %>%
    ggplot(aes(x=value)) %>%
    + geom_histogram(aes(y=..density..),
                     binwidth=.01) %>%
    + facet_grid(~ variable) %>%
    + theme_bw() %>%
    plot

dev.off()

}
