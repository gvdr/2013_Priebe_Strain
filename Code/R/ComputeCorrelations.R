rm(list = ls())
# set ourselves at the root of the Projects folder
setwd("/users/math/gvd16/Projects/2013_Priebe_Strain")
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

Sdf_name <- "Data/Elaborated/Strains_3m.csv"
Mdf_name <- "Data/Elaborated/Mean_distances.csv"
Strain_df <- read.csv(Sdf_name, head = TRUE)
Md_df <- read.csv(Mdf_name, head = TRUE)

correls_diag_corr <- data.frame(Food_Web = character(),
                      Mode = character(),
                      Rank = numeric(),
                      r = numeric(),
                      P = numeric())
for(nam in fws){
  for(ran in 1: Max_rank){
    for(mod in c("in", "out", "all")){
      if((mod != "all") & (ran == 1)) next()
      if(mod == "all") mod_inv <- "all"
      if(mod == "in") mod_inv <- "out"
      if(mod == "out") mod_inv <- "in"
      x_temp <- Strain_df[(Strain_df$Food_Web == nam) & (Strain_df$Rank == ran) & (Strain_df$Mode == mod),"Strain"]
      y_temp <- Md_df[(Md_df$Food_Web == nam) & (Md_df$Rank == ran) & (Md_df$Mode == mod_inv),"Mean_dist"]
      print(paste("Name =", nam, "Rank =", ran, "Mode =", mod, "Dim Rank =", length(x_temp), "Dim meandist =", length(y_temp), sep = " "))
      temp_corr <- cor.test(x_temp,y_temp)
      X <-  data.frame(Food_Web = nam,
                       Mode = mod,
                       Rank = ran,
                       r = temp_corr$estimate,
                       P = temp_corr$p.value)
      correls_diag_corr <- rbind(correls_diag_corr,X)
    }       
  }
}



base_plot_r <- ggplot(correls_diag_corr, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1)
geom_plot_r <- base_plot_r + facet_grid( . ~ Food_Web) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + theme_white_hc
base_plot_P <- ggplot(correls_diag_corr, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1)
geom_plot_P <- base_plot_P + facet_grid( . ~ Food_Web) + scale_color_viridis() #  + guides(colour = "none")
final_plot_P <- geom_plot_P + geom_hline(yintercept=0.05, colour = "red", linetype = 2) + theme_white_hc
all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
all_plot

pdf("Data/Elaborated/MdS_corr_rP.pdf",10,10, useDingbats= FALSE)
all_plot

dev.off()

correls_diag_all <- correls_diag_corr[correls_diag_corr$Mode == "all", ] 
base_plot_r <- ggplot(correls_diag_all, aes(x = Rank)) + geom_line(aes(y = r), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
geom_plot_r <- base_plot_r + facet_grid( . ~ Food_Web) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + theme_white_hc
pdf("Data/Elaborated/MdS_corr_all_r.pdf",10,10, useDingbats= FALSE)
final_plot_r
dev.off()

base_plot_r <- ggplot(correls_diag_all, aes(x = Rank)) + geom_line(aes(y = P), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
geom_plot_r <- base_plot_r + facet_grid( . ~ Food_Web) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + theme_white_hc
pdf("Data/Elaborated/MdS_corr_all_P.pdf",10,10, useDingbats= FALSE)
final_plot_r
dev.off()