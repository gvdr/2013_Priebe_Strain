img_for_corr <- function(obj, col, title, lims=c(0,1)){
   image(obj, col=col, zlim= lims, main = title)
}

correre <- function(Mdf, value, colors, plots=TRUE, mins=TRUE){
  S <- dim(Mdf)[1] / 15
  x <- matrix(nrow = S, ncol = 15)
  for(i in 1:15){
    x[,i] <- Mdf[(Mdf$Rank == i),value]
  }
  Corr <- rcorr(x,type="spearman")
  if(plots == TRUE){
  img_for_corr(Corr$r, colors, paste(Mdf$Food_Web[1],Mdf$Rank[1],sep = " "))
  }
  if(mins == TRUE){
    print(paste("mean = ", mean(Corr$r)))
  }
  return(Corr)
}
 
pdf(file = "Data/Elaborated/Md_corr.pdf",width= 20, height = 10,useDingbats=F)
par(mfrow = c(3,5))
for(mod in c("all","in","out")){
  for(nam in fws){
    Temp_df <- Mdf[(Mdf$Food_Web == nam) & (Mdf$Mode == mod),c("Food_Web","Species","Rank","Mean_dist")]
    n_col <- dim(Temp_df)[1] / 15
    colors <- dichromat(viridis(n_col), "protan")
    correre(Temp_df,"Mean_dist",colors)
  }
}
dev.off()
