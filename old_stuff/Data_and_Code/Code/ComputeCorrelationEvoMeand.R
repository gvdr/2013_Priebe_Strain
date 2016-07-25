rm(list = ls())
options(stringsAsFactors=FALSE)
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
fws_trees <- character()
for(nam in fws){
  if(file.exists(paste("Data/Raw/",nam,"/tree.phy", sep = ""))){
    fws_trees <- cbind(fws_trees,nam)
  }
}


Sdf_name <- "Data/Elaborated/Strains.csv"
Mdf_name <- "Data/Elaborated/Mean_distances.csv"
Ev_name <- "Data/Elaborated/Ev_df.csv"
Strain_df <- read.csv(Sdf_name)
Md_df <- read.csv(Mdf_name)
Ev_df <- read.csv(Ev_name)

correls_meandist <- data.frame(Food_Web = character(),
                               Measure = character(),
                               Rank = numeric(),
                               r = numeric(),
                               P = numeric())

for(nam in fws_trees){
  nam <- "Serengeti_Baskerville"
  tree <- read.newick(paste("Data/Raw/",nam,"/tree.phy",sep = ""))
  tree <- collapse.singles(tree)
  for(meas in c("PE","FP")){
    for(ran in 1: Max_rank){
      temp_s_df <- Md_df[(Md_df$Food_Web == nam) & (Md_df$Rank == ran) & (Md_df$Mode == "all"),]
      temp_ev_df <- Ev_df[(Ev_df$Food_Web == nam) & (Ev_df$Measure == meas),]
      temp_merge <- merge(temp_s_df, temp_ev_df, "Species")
      rownames(temp_merge) <- temp_merge[,"Species"]
      trt <- treedata(tree, temp_merge, sort=T, warnings=F)
      data_temp <- data.frame(Mean_dist = scale(as.numeric(trt$data[,"Mean_dist"])),
                              Distinctiveness = scale(as.numeric(trt$data[,"Distinctiveness"])))
#       data_temp <- data.frame(Mean_dist = rank( - as.numeric(trt$data[,"Mean_dist"]), ties = "average"),
#                               Distinctiveness = rank( - as.numeric(trt$data[,"Distinctiveness"]), ties = "average"))
      rownames(data_temp) <- trt$phy$tip.label
      plot(data_temp, main = paste(meas, ran))
#       temp_corr <- gls(Mean_dist ~ Distinctiveness, data = data_temp, correlation=corPagel(1,trt$phy, fixed = FALSE))
#       temp_corr_0 <- gls(Strain ~ 1, data = data_temp, correlation=corBrownian(1,trt$phy))
#       temp_corr.m1 <- update(temp_corr, . ~ ., method = "ML")
#       temp_corr_0.m1 <- update(temp_corr_0, . ~ ., method = "ML")
#       temp_corr <- phylolm(Strain ~ Distinctiveness, data = data_temp, phy =  trt$phy)
#       temp_corr <- rcorr(data_temp[,"Mean_dist"],data_temp[,"Distinctiveness"], type="spearman")
      temp_corr <- lm(Mean_dist ~ 0 + Distinctiveness, data = data_temp)
      X <-  data.frame(Food_Web = nam,
                       Measure = meas,
                       Rank = ran,
                       
#                        r = temp_corr$r[1,2],
#                        P = temp_corr$P[1,2])
#                        P = anova(temp_corr_0.m1,temp_corr.m1)$p[2])
      correls_meandist <- rbind(correls_meandist,X)
    }
  }
}

base_plot_r <- ggplot(correls_meandist, aes(x = Rank, y = r, colour = Measure)) + geom_line(size = 1)
geom_plot_r <- base_plot_r + facet_grid( . ~ Food_Web) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + theme_white_hc
base_plot_P <- ggplot(correls_meandist, aes(x = Rank, y = P, colour = Measure)) + geom_line(size = 1)
geom_plot_P <- base_plot_P + facet_grid( . ~ Food_Web) + scale_color_viridis() #  + guides(colour = "none")
final_plot_P <- geom_plot_P + theme_white_hc
all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
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