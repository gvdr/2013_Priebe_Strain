##############
rm(list=ls())
# set ourselves at the root of the Projects folder
setwd("/users/math/gvd16/Projects/2013_Priebe_Strain")
options(stringsAsFactors = FALSE)
source("Code/R/LoadPackages.R")
#load the libraries involved in parallel computation
#and start a cluster?
load_packages(parallel = FALSE)
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

###################
adj_t <- Adjs[[3]]
N <- dim(adj_t)[1]
nspecies <- rownames(adj_t)

modes <- c("in", "out", "all")
maxranker <- function(x) ifelse(x %in% c("in", "out"), 7, 4)
L_total <- N * (2 * 7 + 4)
convex_hull_contribution <- data.frame(Food_Web = character(L_total),
                                           Species = character(L_total),
                                           Mode = character(L_total),
                                           Rank = numeric(L_total),
                                           Area_diffs = numeric(L_total),
                                           Vol_diffs = numeric(L_total))

count <- 1
for(mode in modes){
  max_rank <- maxranker(mode)
  space_all <- ASFWE(adj_t,max_rank)
  for(ran in 2:max_rank){
    if(mode == "all"){
      Vars = c(seq(1,ran),seq(max_rank+1, max_rank + ran))
    } else if(mode == "in") {
      Vars = seq(1,ran)
    } else if(mode == "out") {
      Vars = seq(max_rank,max_rank + ran)
    } 
    conv_all <- convhulln(space_all[,Vars], options = "FA")
    area_all <- conv_all$area
    vol_all <- conv_all$vol
    area_diffs <- vol_diffs <- numeric(N)
    start <- count
    for(i in 1:N){
      conv_i <- convhulln(space_all[-i,Vars], options = "FA")
      area_diffs[i] <- area_all - conv_i$area
      vol_diffs[i] <- vol_all - conv_i$vol
      count <- count + 1
      print(paste("mode", mode, "rank", ran))
    }
    end <- count - 1
    convex_hull_contribution[start:end,] <- data.frame(Food_Web = as.character(rep("Serengeti_Baskerville",N)),
                                                           Species = as.character(nspecies),
                                                           Mode = as.character(rep(mode,N)),
                                                           Rank = rep(ran,N),
                                                           Area_diffs = area_diffs,
                                                           Vol_diffs = vol_diffs)
    
  }
}

write.csv(file = "Data/Elaborated/Contributions_to_FD_SB.csv", x = convex_hull_contribution, row.names = FALSE)

rm(list=ls())
# set ourselves at the root of the Projects folder
setwd("/users/math/gvd16/Projects/2013_Priebe_Strain")
options(stringsAsFactors = FALSE)
source("Code/R/LoadPackages.R")
#load the libraries involved in parallel computation
#and start a cluster?
load_packages(parallel = FALSE)
#load the main functions involved in the strain computation
source("Code/R/strain_functions.R")
source("Code/R/ggset.R")


CtFD_name <- "Data/Elaborated/Contributions_to_FD_SB.csv"
Sdf_name <- "Data/Elaborated/Strains_3m.csv"

convex_hull_contribution <- read.csv(CtFD_name, header = TRUE)
Strain_df <- read.csv(Sdf_name, head = TRUE)

correls_strain_FD <- data.frame(Food_Web = character(),
                                Mode = character(),
                                Rank = numeric(),
                                C_Measure = character(),
                                r = numeric(),
                                P = numeric())

for(mod in c("in", "out", "all")){
  Max_rank_mode <- max(convex_hull_contribution[(convex_hull_contribution$Mode == mod),"Rank"])
  for(C_meas in c("Area_diffs", "Vol_diffs")){
    for(ran in 2:Max_rank_mode){
      x_temp <- Strain_df[(Strain_df$Food_Web == "Serengeti_Baskerville") & (Strain_df$Rank == ran) & (Strain_df$Mode == mod),"Strain"]
      y_temp <-  convex_hull_contribution[(convex_hull_contribution$Rank == ran) & (convex_hull_contribution$Mode == mod),C_meas]
      temp_corr <- cor.test(x_temp,y_temp)
      X <-  data.frame(Food_Web = "Serengeti_Baskerville",
                       Mode = mod,
                       Rank = ran,
                       C_Measure = C_meas,
                       r = temp_corr$estimate,
                       P = temp_corr$p.value)
      correls_strain_FD <- rbind(correls_strain_FD,X)
    }
  }       
}

#####
dev.off()
correls_strain_voldiffs <- correls_strain_FD[correls_strain_FD$C_Measure == "Vol_diffs", ]
base_plot_r <- ggplot(correls_strain_voldiffs, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1)
geom_plot_r <- base_plot_r + scale_color_viridis() + coord_cartesian(ylim=c(0, 1))# + guides(colour = "none")
final_plot_r <- geom_plot_r + geom_hline(yintercept=0, colour = "red", linetype = 2) + theme_white_hc

base_plot_P <- ggplot(correls_strain_voldiffs, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1)
geom_plot_P <- base_plot_P + scale_color_viridis() + coord_cartesian(ylim=c(-0.1, 1)) #  + guides(colour = "none")
final_plot_P <- geom_plot_P + geom_hline(yintercept=0.05, colour = "red", linetype = 2) + theme_white_hc

all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
plot(all_plot)

pdf("Data/Plots/SB_Strain_FD.pdf",2,4, useDingbats= FALSE)
plot(all_plot)
dev.off()
######

Mdf_name <- "Data/Elaborated/Mean_distances.csv"
Md_df <- read.csv(Mdf_name, head = TRUE)


#####
correls_meands_FD <- data.frame(Food_Web = character(),
                                Mode = character(),
                                Rank = numeric(),
                                C_Measure = character(),
                                r = numeric(),
                                P = numeric())

for(mod in c("in", "out", "all")){
  if(mod == "all") mod_inv <- "all"
  if(mod == "in") mod_inv <- "out"
  if(mod == "out") mod_inv <- "in"
  Max_rank_mode <- max(convex_hull_contribution[(convex_hull_contribution$Mode == mod),"Rank"])
  for(C_meas in c("Area_diffs", "Vol_diffs")){
    for(ran in 2:Max_rank_mode){
      x_temp <- Md_df[(Md_df$Food_Web == "Serengeti_Baskerville") & (Md_df$Rank == ran) & (Md_df$Mode == mod_inv),"Mean_dist"]
      y_temp <-  convex_hull_contribution[(convex_hull_contribution$Rank == ran) & (convex_hull_contribution$Mode == mod),C_meas]
      temp_corr <- cor.test(x_temp,y_temp)
      X <-  data.frame(Food_Web = "Serengeti_Baskerville",
                       Mode = mod,
                       Rank = ran,
                       C_Measure = C_meas,
                       r = temp_corr$estimate,
                       P = temp_corr$p.value)
      correls_meands_FD <- rbind(correls_meands_FD,X)
    }
  }       
}

correls_meands_voldiffs <- correls_meands_FD[correls_meands_FD$C_Measure == "Vol_diffs", ]
base_plot_r_m <- ggplot(correls_meands_voldiffs, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1)
geom_plot_r_m <- base_plot_r_m + scale_color_viridis() + coord_cartesian(ylim=c(0, 1)) # + guides(colour = "none")
final_plot_r_m <- geom_plot_r_m + geom_hline(yintercept=0, colour = "red", linetype = 2) + theme_white_hc

base_plot_P_m <- ggplot(correls_meands_voldiffs, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1)
geom_plot_P_m <- base_plot_P_m + scale_color_viridis()  + coord_cartesian(ylim=c(-0.1, 1))#  + guides(colour = "none")
final_plot_P_m <- geom_plot_P_m + geom_hline(yintercept=0.05, colour = "red", linetype = 2) + theme_white_hc

all_plot_m <- grid_arrange_shared_legend(final_plot_r_m,final_plot_P_m)
plot(all_plot_m)

pdf("Data/Plots/SB_Meands_FD.pdf",2,4, useDingbats= FALSE)
plot(all_plot_m)
dev.off()