rm(list = ls())
# set ourselves at the root of the Projects folder
setwd("/users/math/gvd16/Projects/2013_Priebe_Strain")
source("Code/R/LoadPackages.R")
#load the libraries involved in parallel computation
#and start a cluster?
load_packages(parallel = FALSE)
#load the function we need
source("Code/R/ed_arne_functions.R")
#retrieve all the food webs
fws <- list.dirs("Data/Elaborated/interactions")
fws_trees <- character()
for(nam in fws){
  if(file.exists(paste("Data/Raw/",nam,"/tree.phy", sep = ""))){
    fws_trees <- cbind(fws_trees,nam)
  }
}

Ev_df <- data.frame(Food_Web = character(),
                    Measure = character(),
                    Species = character(),
                    Distinctiveness = numeric(),
                    stringsAsFactors=F)
#Ev_list <- vector('list', 3)
#get the trees
for(i in 1:3){
  nam <- fws_trees[i]
  tree_path <- paste("Data/Raw/",nam,"/tree.phy", sep = "")
  fp_path <- paste("Data/Elaborated/evo_dist/",nam,"_fp.csv", sep = "")
  pe_path <- paste("Data/Elaborated/evo_dist/",nam,"_pe.csv", sep = "")
  Ev_temp <- compute_fp_ep(tree_path,fp_path,pe_path)
  N <- dim(Ev_temp$Equal_Split)[1]
  Species_nam <- as.character(rownames(Ev_temp$Equal_Split))
  X <- data.frame(Food_Web = rep(as.character(nam),N),
                  Measure = rep("PE",N),
                  Species = Species_nam,
                  Distinctiveness = Ev_temp$Equal_Split[,1],
                  stringsAsFactors=F)
  Y <- data.frame(Food_Web = rep(as.character(nam),N),
                  Measure = rep("FP",N),
                  Species = Species_nam,
                  Distinctiveness = Ev_temp$Fair_Proportion[,1],
                  stringsAsFactors=F)
  rownames(X) <- rownames(Y) <- NULL
  Ev_df <- rbind(Ev_df,rbind(X,Y))
  #Ev_list[[i]] <- Ev_temp
}

write.csv(file = "Data/Elaborated/Ev_df.csv", x = Ev_df, row.names = FALSE)
