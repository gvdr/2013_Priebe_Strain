####### Effect of energy flows on RDPG

nb_cores <- 12L
Mean_Noise <- 0.5
Var_Noise <- 0.5
NSIM <- 100

RANK <- 3

###############
setwd("/users/math/gvd16/Projects/2013_Priebe_Strain/WD")
options(stringsAsFactor = FALSE)
source("../Code/R/LoadPackages.R")
load_packages(parallel = TRUE)
source("../Code/R/PlotFunctions_for_webs.R")
source("../Code/R/ggset.R")
source("../Code/R/Noiser.R")
library(broom)
library(parallel)
library(grid)
#################
#retrieve all the food webs
fws <- list.dirs("../Data/Elaborated/interactions")
M <- length(fws)

Adjs <- vector('list',M)
for(i in 1:M){
  #retrieve all adjacency matrix
  fws[i] %>% get_Adj(., base = paste("../Data/Elaborated/interactions/",.,sep="")) -> Adjs[[i]]
}
###############
### Loop through nets and modes and plot
for(i in 1:5){
  TheChosenNetwork <- fws[i]
  adj <- Adjs[[i]]
  for(MODE in c("inward","outward","total")){
    
    if(MODE == "inward") asfwe <- ASFWE_in
    if(MODE == "outward") asfwe <- ASFWE_out
    if(MODE == "total") asfwe <- ASFWE
    
    adj %>% asfwe(RANK) %>%
      dist %>% as.matrix %>%
      rowMeans %>% setNames(row.names(adj))-> topo_meand

    
    mclapply( 1:NSIM,
              function(x){
                comparison_flow_meand(adj,Rank=RANK,Mode=MODE,Mean_Noise,Var_Noise)
              },
              mc.cores =  getOption("mc.cores", nb_cores)
    ) -> ot
    
    sapply(1:NSIM,
           function(x){
             ot[[x]] %>%
               .[names(topo_meand)] %>%
{lm(.~ 0 + topo_meand)} -> model
c(model %>% tidy %$% p.value,
  model %>% summary %$% adj.r.squared)
           } %>%
  setNames(c("p.value","adj.r.squared"))
    ) %>% t -> simsvalues

How_Many <- 4
Higher_Window <- 10

topo_meand %>%
  sort(decreasing = TRUE) %>%
  names %>%
  .[1:How_Many] -> topo_higher

sapply(1:NSIM,
       function(x){
         ot[[x]] %>%
           sort(decreasing = TRUE) %>%
           names %>%
           .[1:Higher_Window] -> higher_strain
         as.numeric(topo_higher %in% higher_strain) %>%
           setNames(topo_higher) 
       }
) %>% t -> freq_higher

colnames(freq_higher) %<>% gsub(x=.,"_"," ")

freq_higher %>% colMeans %>%
  as.data.frame %>%
  setNames("Frequency") %>%
  mutate(species_names = colnames(freq_higher)) %>%
  ggplot(aes(x=reorder(species_names, -topo_meand[species_names]),y=Frequency)) +
  geom_point() +
  ylim(c(0,1)) +
  ggtitle(TheChosenNetwork %>% gsub(x=.,"_"," ") %>% paste0(" ", MODE)) +
  labs(x="Taxa id",
       y=paste0("Frequency in top ",How_Many)) +
  geom_hline(yintercept = 0.5,color="green",linetype = 2) +
  geom_hline(yintercept = Higher_Window/length(topo_meand),color="red",linetype = 3) +
  theme_minimal() +
  theme(axis.title.x = element_text(vjust=9),
        axis.text.x = element_text(size=10,angle=45, vjust=4,face = 'italic'),
        plot.margin = unit(c(1,1,-2,1),"cm") ) -> p

paste0(TheChosenNetwork,"_",MODE,".pdf") %>%
  ggsave(p,path="../Plots/Weighted/meand/")



simsvalues %>% .[,1] %>% "<"(0.05) %>% sum %>% "/"(NSIM) -> freq.p
simsvalues %>% .[,2] %>% mean -> meanr
paste0(TheChosenNetwork," ",MODE,", significants = ", freq.p %>% percent,
       ", mean r squared= ", meanr %>% percent) %>%
  print
paste0("../Plots/Weighted/meand/",TheChosenNetwork,"_bp_",MODE,".pdf") %>%
  pdf(width=6,height=6,useDingbats = FALSE)
simsvalues %>% boxplot
points(c(0.05,0.5),col="red",pch="_",cex=5)
dev.off()
  }
}
