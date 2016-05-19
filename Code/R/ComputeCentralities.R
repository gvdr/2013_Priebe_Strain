fws <- list.dirs("Data/Elaborated/interactions")
M <- length(fws)
#Initialize a list for the adjacency matrices
Gs <- vector('list',M)
for(i in 1:M){
  #retrieve all adjacency matrix
  fws[i] %>% get_Graph(., base = paste("Data/Elaborated/interactions/",.,sep="")) -> Gs[[i]]
}


Centralities <- data.frame(Food_Web = character(),
                           Species = character(),
                           DC = numeric(),
                           BC = numeric(),
                           CC = numeric(),
                           EC = numeric(),
                           IC = numeric(),
                           SC = numeric())

for(i in 1:M){
  name_temp <- fws[i]
  graph_temp <- Gs[[i]]
  cent_temp <- compute_centralities(graph_temp,name_temp)
  Centralities <- rbind(Centralities,cent_temp)
}