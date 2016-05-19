mod <- "all"
min_rank <- 1 #1 for "all", 2 for "in" and "out"
max_rank <- 15

mat_dim <- max_rank - (min_rank - 1)
Correls_strain_P <- matrix(NA, mat_dim, mat_dim)
Correls_strain_r <- matrix(NA, mat_dim, mat_dim)
#, "out", "all")

for(rani in min_rank:max_rank){
    straini <- Strain_df[(Strain_df$Food_Web == "Serengeti_Baskerville") & (Strain_df$Mode == mod) & (Strain_df$Rank == rani), "Strain"]
    straini <- rank(straini, tie = "average")
    ii <- rani - (min_rank - 1)
    for(ranj in min_rank:max_rank){
      strainj <- Strain_df[(Strain_df$Food_Web == "Serengeti_Baskerville") & (Strain_df$Mode == mod) & (Strain_df$Rank == ranj), "Strain"]
      strainj <- rank(strainj, tie = "average")
      jj <- ranj - (min_rank - 1)
      Correls_strain_P[ii,jj] <- rcorr(straini,strainj)$P[1,2]
      Correls_strain_r[ii,jj] <- rcorr(straini,strainj)$r[1,2]
    }
  }

Data_corr_P <- melt(Correls_strain_P) 
Data_corr_P$Var1 <- Data_corr_P$Var1 + (min_rank - 1)
Data_corr_P$Var2 <- Data_corr_P$Var2 + (min_rank - 1)
Data_corr_P$Mode <- mod
Data_corr_P$stat <- "P"
names(Data_corr_P) <- c("Rank_x", "Rank_y", "value", "Mode", "Statistics")
heat_plot_P <- ggplot(Data_corr_P, aes(factor(Rank_x), factor(Rank_y))) + geom_tile(aes(fill = value)) + labs(title= paste(mod, "P"), x = "Rank", y = "Rank")
heat_plot_P_colors <- heat_plot_P + theme_white_hc  + scale_fill_viridis(discrete = FALSE, trans = "reverse", limits = c(0.05, 0.0))#, limits=c(-10, 10), breaks=seq(0, 0.05, by=0.01))

Data_corr_r <- melt(Correls_strain_r) 
Data_corr_r$Var1 <- Data_corr_r$Var1 + (min_rank - 1)
Data_corr_r$Var2 <- Data_corr_r$Var2 + (min_rank - 1)
Data_corr_r$Mode <- mod
Data_corr_r$stat <- "P"
names(Data_corr_r) <- c("Rank_x", "Rank_y", "value", "Mode", "Statistics")
heat_plot_r <- ggplot(Data_corr_r, aes(factor(Rank_x), factor(Rank_y))) + geom_tile(aes(fill = value)) + labs(title= paste(mod, "r"), x = "Rank", y = "Rank")
heat_plot_r_colors <- heat_plot_r + theme_white_hc  + scale_fill_viridis(discrete = FALSE, limits = c(0.0, 1.0))

plottone <- grid.arrange(heat_plot_P_colors,heat_plot_r_colors, ncol=2)


rani <- 3
straini <- Strain_df[(Strain_df$Food_Web == "Serengeti_Baskerville") & (Strain_df$Mode == mod) & (Strain_df$Rank == rani), c("Species","Strain")]

straini_sorted <- sort.data.frame(straini, by = "Strain")
rownames(straini_sorted) <- seq(length=nrow(straini_sorted))
straini_sorted$Species <- factor(straini_sorted$Species, levels = straini_sorted$Species)

plot_s <- ggplot(data = straini_sorted) +
  geom_point(aes(y = Strain, x = Species)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = .5, color = "lightgrey"),
        panel.grid.major = element_blank(),
        axis.line = element_line(size=.7, color = "grey"),
        text = element_text(size=14),
        strip.background = element_rect(fill= NA, color = NA),
        strip.text = element_text(size=14, family = "Helvetica"),
        axis.text.x=element_blank()) +
  scale_y_log10()
plot_s
