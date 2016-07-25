rm(list=ls())
# set ourselves at the root of the Projects folder
setwd("/users/math/gvd16/Projects/2013_Priebe_Strain")
options(stringsAsFactor = FALSE)
source("Code/R/LoadPackages.R")
#load the libraries involved in parallel computation
#and start a cluster?
load_packages(parallel = TRUE)
#load the main functions involved in the strain computation
source("Code/R/strain_functions.R")
source("Code/R/ggset.R")
#set the Maximum rank we are going to analyse
Max_rank <- 15 
#retrieve all the food webs
fws <- list.dirs("Data/Elaborated/interactions")
M <- length(fws)
#Initialize a list for the adjacency matrices
Adjs <- vector('list',M)
for(i in 1:M){
  #retrieve all adjacency matrix
  fws[i] %>% get_Adj(., base = paste("Data/Elaborated/interactions/",.,sep="")) -> Adjs[[i]]
}

Strain_df <- data.frame(
			Food_Web = character(),
      Mode = character(),
			Species = character(),
			Rank = numeric(),
			Strain = numeric()
			)

for(name in list.dirs("Data/Elaborated/interactions")){
  name %>%
    get_Adj( . , base = paste("Data/Elaborated/interactions/",.,sep="")) -> Adj
  N <- dim(Adj)[1]
  labels <- rownames(Adj)
  fw <- rep(name,N)
  for(rank in 1:max_rank){
    for(mod in c("out", "in", "all")){
      #compute strain
      rank_index <- rep(rank,N)
      mods <- rep(mod,N)
      strains <- compute_strain(Adj,rank,mod)
      Strain_temp <- data.frame(Food_Web = fw,
                                Mode = mods,
                                Species = labels,
                                Rank = rank_index,
                                Strain = strains)
      Strain_df <- rbind(Strain_df,Strain_temp)
    }
  }
}

write.csv(file = "Data/Elaborated/Strains_3m.csv", x = Strain_df, row.names = FALSE)

#save to an Rdata object for data analysis and visualization
save(Strain_df, file="Elaborated/Strain_by_rank/sAll.Rdata")

base_plot <- ggplot(Strain_df, aes(x = Rank, y = Strain, colour = factor(Species))) + geom_line(size = 0.4) + scale_color_manual(values = sample(viridis(933), replace = FALSE))
geom_plot <- base_plot + facet_grid(Mode ~ Food_Web, scales = "free") + guides(colour = "none") #values=rep(brewer.pal(5,"Set2"),times=211))
themed_plot <- geom_plot + theme_white_hc + scale_y_log10(limits = c(10^(-8),NA))
themed_plot
ggsave(themed_plot, file="strains_log_cut.pdf", height=10, width=20, units='in', dpi=600)
