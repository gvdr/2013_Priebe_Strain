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
source("Code/R/AiccPhylo.R")
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

###Load DataSets
#######
Sdf_name <- "Data/Elaborated/Strains_3m.csv"
Mdf_name <- "Data/Elaborated/Mean_distances.csv"
Ev_name <- "Data/Elaborated/Ev_df.csv"
Strain_df <- read.csv(Sdf_name, colClasses=c("Food_Web" = "character",
                                             "Mode" = "character",
                                             "Species" = "character",
                                             "Rank" = "numeric",
                                             "Strain" = "numeric"))
Md_df <- read.csv(Mdf_name, colClasses=c("Food_Web" = "character",
                                         "Mode" = "character",
                                         "Rank" = "numeric",
                                         "Species" = "character",
                                         "Mean_dist" = "numeric"))
Ev_df <- read.csv(Ev_name, colClasses=c("Food_Web" = "character",
                                        "Measure" = "character",
                                        "Species" = "character",
                                        "Distinctiveness" = "numeric"))
##########################
Correls_MS <- data.frame(Food_Web = character(),
                         Mode = character(),
                         Rank = numeric(),
                         Sigma = numeric(),
                         P = numeric())

Meandist_aiccs <- data.frame(Food_Web = character(),
                             Mode = character(),
                             Rank = numeric(),
                             ZeroAICC = numeric(),
                             BmAICC = numeric(),
                             OuAICC = numeric(),
                             EbAICC = numeric())

Strain_aiccs <- data.frame(Food_Web = character(),
                           Mode = character(),
                           Rank = numeric(),
                           ZeroAICC = numeric(),
                           BmAICC = numeric(),
                           OuAICC = numeric(),
                           EbAICC = numeric())
###########################################
for(nam in c("Serengeti_deVisser","Serengeti_Baskerville")){
  tree_temp <- read.newick(paste("Data/Raw/",nam,"/tree.phy",sep = ""))
  tree_temp <- collapse.singles(tree_temp)
  tree_temp <- multi2di(tree_temp)
  for(mod in c("out","in", "all")){
      for(ran in 1: Max_rank){
        if((mod != "all") & (ran == 1)) next()
        if(mod == "all") mod_inv <- "all"
        if(mod == "in") mod_inv <- "out"
        if(mod == "out") mod_inv <- "in"
        temp_md_df <- Md_df[(Md_df$Food_Web == nam) &
                              (Md_df$Rank == ran) &
                              (Md_df$Mode == mod),
                            c("Species","Mean_dist")]
        
        temp_s_df <- Strain_df[(Strain_df$Food_Web == nam) &
                                 (Strain_df$Rank == ran) &
                                 (Strain_df$Mode == mod_inv),
                               c("Species","Strain")]
        
        temp_merge <- merge(temp_md_df, temp_s_df,"Species")
        rownames(temp_merge) <- temp_merge[,"Species"]
        temp_merge <- temp_merge[,c("Strain","Mean_dist")]
        #match phylogeny and dataset
        trt <- treedata(tree_temp, temp_merge, warnings=F)
        fitStrain <- AiccWeights(trt,"Strain")
        fitMeanDist <- AiccWeights(trt,"Mean_dist")
        
        pglsModel <- gls(Strain ~ Mean_dist,
                         data = as.data.frame(trt$data),
                         correlation = corMartins((fitStrain$alpha + fitMeanDist$alpha)/2,
                                                  phy = trt$phy,
                                                  fixed = T),
                         method = "ML")
        
        Co <- data.frame(Food_Web = nam,
                         Mode = mod,
                         Rank = ran,
                         Sigma = summary(pglsModel)$sigma,
                         P = anova(pglsModel)$p[2])
        
#         Sa <- data.frame(Food_Web = nam,
#                          Mode = mod,
#                          Rank = ran,
#                          ZeroAICC = fitStrain$weights[1],
#                          BmAICC = fitStrain$weights[2],
#                          OuAICC = fitStrain$weights[3],
#                          EbAICC = fitStrain$weights[4])
# 
#         Ma <- data.frame(Food_Web = nam,
#                          Mode = mod,
#                          Rank = ran,
#                          ZeroAICC = fitMeanDist$weights[1],
#                          BmAICC = fitMeanDist$weights[2],
#                          OuAICC = fitMeanDist$weights[3],
#                          EbAICC = fitMeanDist$weights[4])

        Correls_MS <- rbind(Correls_MS,Co)
#         Strain_aiccs <- rbind(Strain_aiccs,Sa)
#         Meandist_aiccs <- rbind(Meandist_aiccs,Ma)
      }
    }
  }
###########################################
base_plot_r <- ggplot(Correls_MS, aes(x = Rank, y = r, colour = Mode)) + geom_line(size = 1) + geom_point(size = 3)
geom_plot_r <- base_plot_r + facet_grid( . ~ Food_Web) + scale_color_viridis() + geom_hline(yintercept=0) # + guides(colour = "none")
final_plot_r <- geom_plot_r + theme_white_hc
base_plot_P <- ggplot(Correls_MS, aes(x = Rank, y = P, colour = Mode)) + geom_line(size = 1) + geom_point(size = 3)
geom_plot_P <- base_plot_P + facet_grid( . ~ Food_Web) + scale_color_viridis() + geom_hline(yintercept=0.05) #  + guides(colour = "none")
final_plot_P <- geom_plot_P + theme_white_hc
all_plot <- grid_arrange_shared_legend(final_plot_r,final_plot_P)
pdf("Data/Elaborated/MdS_corr_P.pdf",10,10, useDingbats= FALSE)
all_plot
dev.off()


base_plot_m <- ggplot(Meandist_aiccs_melted[Meandist_aiccs_melted$Food_Web == "Serengeti_Baskerville", ], aes(variable, value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = as.factor(Rank)),position = position_jitter(width = .3), size = 5)#aes(y = r), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
geom_plot_m <- base_plot_m + facet_grid( . ~ Mode) + scale_color_viridis(alpha = 0.8) # + guides(colour = "none")
final_plot_m <- geom_plot_m + theme_white_hc
final_plot_m
base_plot_s <- ggplot(Strain_aiccs_melted[Strain_aiccs_melted$Food_Web == "Serengeti_Baskerville", ], aes(variable, value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = as.factor(Rank)),position = position_jitter(width = .3), size = 5)#aes(y = r), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
geom_plot_s <- base_plot_s + facet_grid( . ~ Mode) + scale_color_viridis(alpha = 0.8) # + guides(colour = "none")
final_plot_s <- geom_plot_s + theme_white_hc
final_plot_s


abase_plot_m <- ggplot(Meandist_aiccs_melted, aes(variable, value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = as.factor(Rank)),position = position_jitter(width = .3), size = 5)#aes(y = r), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
ageom_plot_m <- abase_plot_m + facet_grid( Food_Web ~ Mode) + scale_color_viridis(alpha = 0.8) # + guides(colour = "none")
afinal_plot_m <- ageom_plot_m + theme_white_hc
afinal_plot_m
abase_plot_s <- ggplot(Strain_aiccs_melted, aes(variable, value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = as.factor(Rank)),position = position_jitter(width = .3), size = 5)#aes(y = r), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
ageom_plot_s <- abase_plot_s + facet_grid( Food_Web ~ Mode) + scale_color_viridis(alpha = 0.8) # + guides(colour = "none")
afinal_plot_s <- ageom_plot_s + theme_white_hc
afinal_plot_s

base_plot_all <- ggplot(All_aiccs[All_aiccs$Food_Web == "Serengeti_Baskerville", ], aes(variable, value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = as.factor(Rank)),position = position_jitter(width = .3), size = 5)#aes(y = r), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
geom_plot_all <- base_plot_all + facet_grid( Kind ~ Mode) + scale_color_viridis(alpha = 0.8) # + guides(colour = "none")
final_plot_all <- geom_plot_all + theme_white_hc
final_plot_all

ExportPlot(final_plot_all, "Data/Plots/All_Sig_in", 10,10)

base_plot_r <- ggplot(correls_diag_all, aes(x = Rank)) + geom_line(aes(y = P), size = 1) #+ geom_line(aes(y = P, colour = "p value"), size = 1)
geom_plot_r <- base_plot_r + facet_grid( . ~ Food_Web) + scale_color_viridis() # + guides(colour = "none")
final_plot_r <- geom_plot_r + theme_white_hc
pdf("Data/Elaborated/MdS_corr_all_P.pdf",10,10, useDingbats= FALSE)
final_plot_r
dev.off()

####################

temp_md_df <- Md_df[(Md_df$Food_Web == "Serengeti_Baskerville") &
                      (Md_df$Rank == 3) &
                      (Md_df$Mode == "all"),
                    c("Species","Mean_dist")]

temp_s_df <- Strain_df[(Strain_df$Food_Web == "Serengeti_Baskerville") &
                         (Strain_df$Rank == 3) &
                         (Strain_df$Mode == "all"),
                       c("Species","Strain")]

temp_merge <- merge(temp_md_df, temp_s_df,"Species")
rownames(temp_merge) <- temp_merge[,"Species"]
temp_merge <- temp_merge[,c("Strain","Mean_dist")]

temp_merge$Strain_scaled <- scale(log(temp_merge$Strain))
temp_merge$Mean_dist_scaled <- scale(log(temp_merge$Mean_dist))

trt <- treedata(tree_temp, temp_merge, warnings=F)

Mplot <- contMap(trt$phy,trt$data[,"Mean_dist_scaled"])
Splot <- contMap(trt$phy,trt$data[,"Strain_scaled"])
nM <-length(Mplot$cols)
nS <-length(Splot$cols)
Mplot$cols[ ] <- rev(viridis(nM))
Splot$cols[ ]<- rev(viridis(nS))

CairoSVG(file = "Data/Plots/SB_Strain_vs_Meandist_3.svg")
layout(matrix(1:2,1,2),widths=c(0.5,0.5))
par(cex=1)
plot(Mplot,ftype="off",sig=0,legend=0, outline = F, lwd = 1.5)
plot(Splot,ftype="off",direction="leftwards", sig=0,legend=0, outline = F, lwd = 1.5)
dev.off()