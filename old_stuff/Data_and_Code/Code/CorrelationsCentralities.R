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
C_name <- "Data/Elaborated/Centralities.csv"
Mdf_name <- "Data/Elaborated/Mean_distances.csv"
Strain_df <- read.csv(Sdf_name, head = TRUE)
C_df <- read.csv(C_name, head = TRUE)
Md_df <- read.csv(Mdf_name, head = TRUE)

correls_strain_cs <- data.frame(Food_Web = character(),
                                Mode = character(),
                                Rank = numeric(),
                                Centrality_Measure = character(),
                                r = numeric(),
                                P = numeric())
for(nam in fws){
  for(ran in 1: Max_rank){
    for(mod in c("in", "out", "all")){
      if((mod != "all") & (ran == 1)) next()
      for(meas in c("DC", "BC", "CC", "EC", "IC", "SC")){
        x_temp <- Strain_df[(Strain_df$Food_Web == nam) & (Strain_df$Rank == ran) & (Strain_df$Mode == mod),"Strain"]
        y_temp <- C_df[(C_df$Food_Web == nam) & (C_df$Measure == meas),"Centrality"]
        print(paste("Name =", nam, "Rank =", ran, "Mode =", mod, "Meas =", meas, "Dim Rank =", length(x_temp), "Dim meandist =", length(y_temp), sep = " "))
        temp_corr <- cor.test(x_temp,y_temp)
        X <-  data.frame(Food_Web = nam,
                         Mode = mod,
                         Rank = ran,
                         Centrality_Measure = meas,
                         r = temp_corr$estimate,
                         P = temp_corr$p.value)
        correls_strain_cs <- rbind(correls_strain_cs,X)
      }
    }       
  }
}


correls_strain_sb <- correls_strain_cs[correls_strain_cs$Food_Web == "Serengeti_Baskerville",]
base_plot_r <- ggplot(correls_strain_sb, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1)
geom_plot_r <- base_plot_r + facet_grid( Centrality_Measure ~ . ) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + geom_hline(yintercept=0, colour = "red", linetype = 2) + theme_white_hc
base_plot_P <- ggplot(correls_strain_sb, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1)
geom_plot_P <- base_plot_P + facet_grid( Centrality_Measure ~ .) + scale_color_viridis() #  + guides(colour = "none")
final_plot_P <- geom_plot_P + geom_hline(yintercept=0.05, colour = "red", linetype = 2) + theme_white_hc
all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
all_plot

pdf("Data/Elaborated/SB_Strain_corr_rP.pdf",3,10, useDingbats= FALSE)
all_plot
dev.off()


correls_meands_cs <- data.frame(Food_Web = character(),
                                Mode = character(),
                                Rank = numeric(),
                                Centrality_Measure = character(),
                                r = numeric(),
                                P = numeric())
for(nam in fws){
  for(ran in 1: Max_rank){
    for(mod in c("in", "out", "all")){
      if((mod != "all") & (ran == 1)) next()
      if(mod == "all") mod_inv <- "all"
      if(mod == "in") mod_inv <- "out"
      if(mod == "out") mod_inv <- "in"
      for(meas in c("DC", "BC", "CC", "EC", "IC", "SC")){
        x_temp <- Md_df[(Md_df$Food_Web == nam) & (Md_df$Rank == ran) & (Md_df$Mode == mod_inv),"Mean_dist"]
        y_temp <- C_df[(C_df$Food_Web == nam) & (C_df$Measure == meas),"Centrality"]
        print(paste("Name =", nam, "Rank =", ran, "Mode =", mod, "Meas =", meas, "Dim Rank =", length(x_temp), "Dim meandist =", length(y_temp), sep = " "))
        temp_corr <- cor.test(x_temp,y_temp)
        X <-  data.frame(Food_Web = nam,
                         Mode = mod,
                         Rank = ran,
                         Centrality_Measure = meas,
                         r = temp_corr$estimate,
                         P = temp_corr$p.value)
        correls_meands_cs <- rbind(correls_meands_cs,X)
      }
    }       
  }
}

base_plot_r <- ggplot(correls_meands_cs, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1)
geom_plot_r <- base_plot_r + facet_grid( Centrality_Measure ~ Food_Web ) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + geom_hline(yintercept=0, colour = "red", linetype = 2) + theme_white_hc
base_plot_P <- ggplot(correls_meands_cs, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1)
geom_plot_P <- base_plot_P + facet_grid( Centrality_Measure ~ Food_Web ) + scale_color_viridis() #  + guides(colour = "none")
final_plot_P <- geom_plot_P + geom_hline(yintercept=0.05, colour = "red", linetype = 2) + theme_white_hc
all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
all_plot

pdf("Data/Elaborated/MdS_corr_rP.pdf",10,10, useDingbats= FALSE)
all_plot
dev.off()

correls_meands_cs_sb <- correls_meands_cs[correls_meands_cs$Food_Web == "Serengeti_Baskerville",]
base_plot_r <- ggplot(correls_meands_cs_sb, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1) + coord_cartesian(ylim=c(-0.3, 1))
geom_plot_r <- base_plot_r + facet_grid( Centrality_Measure ~ . ) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + geom_hline(yintercept=0, colour = "red", linetype = 2) + theme_white_hc
base_plot_P <- ggplot(correls_meands_cs_sb, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1)
geom_plot_P <- base_plot_P + facet_grid( Centrality_Measure ~ . ) + scale_color_viridis() #  + guides(colour = "none")
final_plot_P <- geom_plot_P + geom_hline(yintercept=0.05, colour = "red", linetype = 2) + theme_white_hc
all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
all_plot

pdf("Data/Elaborated/SB_MdS_corr_rP.pdf",3,10, useDingbats= FALSE)
all_plot
dev.off()