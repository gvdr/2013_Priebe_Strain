
strain_graph_plot <- function(Name,Graph){
  
  ######################### GRAPH STUFF  #########################
  
  #network nodes statistics
  outd_norm <- normalize.vector(degree(Graph, mode="out"))
  ind_norm <- normalize.vector(degree(Graph, mode="in"))+1
  trophic_level <- (outd_norm)/(ind_norm)
  degree <- degree(Graph, mode="all")
  
  #plot the graph as it is
  graph_title = paste("Graph",Name,sep="_")
  pdf(paste(graph_title,".pdf",sep=""))
  par(mfrow=c(1,1))
  V(Graph)$size <- trophic_level*70
  colbar <- heat.colors(max(degree)+1)
  V(Graph)$color <- colbar[degree]
  plot(Graph, layout=layout.grid, edge.curved=TRUE, edge.width=0.01, vertex.label=NA, edge.arrow.size=0.5)
  dev.off()
  
}


#main plotting magik function
strain_plotting <- function(Name,Graph,Rank){
  
  ######################### GRAPH STUFF  #########################
  
  #basic network analysis
  Size = length(V(Graph))
  #network nodes statistics
  #pr_norm <- normalize.vector(page.rank(Graph)$vector)
  closeness_norm <- normalize.vector(closeness(Graph))
  betweenness_norm <- normalize.vector(betweenness(Graph))
  outd_norm <- normalize.vector(degree(Graph, mode="out"))
  ind_norm <- normalize.vector(degree(Graph, mode="in"))
  alld_norm <- normalize.vector(degree(Graph, mode="all"))
  evcent_norm <- normalize.vector(evcent(Graph,directed=TRUE)$vector)
  
  #########################  ASGE STUFF   #########################
  
  #define the function computing the parameter matrix
  ASFWE = function(A,dim)
  {
    S = svd(A)
    Xhat1 = S$u[,1:dim] %*% diag(sqrt(S$d[1:dim]))
    Xhat2 = S$v[,1:dim] %*% diag(sqrt(S$d[1:dim]))
    Xhat = cbind(Xhat1, Xhat2)
  }
  Adj = get.adjacency(Graph)
  
  #pdf plot the svds
  svd_title = paste("Svds",Name,sep="_")
  pdf(paste(svd_title,".pdf",sep=""))
  par(mfrow=c(1,1))
  plot(svd(Adj)$d^2/sum(svd(Adj)$d^2), log = "x", ylab="Singular Value")
  dev.off()
  
  #compute the rank Rank paramater matrix for the original fw 
  Xhat = ASFWE(Adj,Rank)
  #now calculate the strain (via brute force)
  Vars = seq(1,2*Rank)
  strain = NULL
  for(i in 1:Size)
  {
    Xhat_i = ASFWE(Adj[-i,-i],Rank)
    strain[i] = procOPA( Xhat_i[,Vars] , Xhat[-i,Vars] ,scale=T,reflect=T)$OSS
  }
  
  #pdf pairsplot
  pairs_title = paste("Pairs",Name,Rank,sep="_")
  pdf(paste(pairs_title,".pdf",sep=""))
  par(mfrow=c(1,1))
  pairs(Xhat[,Vars],cex=strain)
  dev.off()
  
  #"asfwe-pairs-influence-families.pdf"
  strain_norm = normalize.vector(strain)
  max_strain_norm = max(strain_norm)
  
  #########################  MAIN PLOTTING  #########################
  
  #define the colors
  color = NULL
  for(i in 1:Size)
  {
    color[i] = 100 * (strain_norm[i]*max_strain_norm) / (max_strain_norm^2)
  }
  colbar <- heat.colors(max(color)+1)
  V(Graph)$color <- colbar[color+1]
  #open the file
  main_title = paste("Plot/Strain_statistics",Name,Rank,sep="_")
  pdf(paste(main_title,".pdf",sep=""), width = 6, height = 8)
  par(mfrow=c(3,2))
  #degree
  V(Graph)$size <- alld_norm^2*400
  plot(Graph, layout=layout.grid, vertex.label=NA, edge.lty=0, edge.arrow.mode=0)
  title(main="degree")
  #outdegree
  V(Graph)$size <- outd_norm^2*400
  plot(Graph, layout=layout.grid, vertex.label=NA, edge.lty=0, edge.arrow.mode=0)
  title(main="outdegree")
  #betweenness
  V(Graph)$size <- betweenness_norm^2*400
  plot(Graph, layout=layout.grid, vertex.label=NA, edge.lty=0, edge.arrow.mode=0)
  title(main="betweenness")
  #closeness
  V(Graph)$size <- closeness_norm^2*400
  plot(Graph, layout=layout.grid, vertex.label=NA, edge.lty=0, edge.arrow.mode=0)
  title(main="closeness")
  #page.rank
  #V(Graph)$size <- pr_norm^2*400
  #plot(Graph, layout=layout.grid, vertex.label=NA, edge.lty=0, edge.arrow.mode=0)
  #title(main="pageranks")
  #eigenvector centrality
  V(Graph)$size <- evcent_norm^2*400
  plot(Graph, layout=layout.grid, vertex.label=NA, edge.lty=0, edge.arrow.mode=0)
  title(main="Eigenvector Centrality")
  #devoff
  dev.off()
  
  #########################  STRAIN DISTRIBUTION #########################
  
  
  
  #pdf plot a histogram of the strain
  strain_title = paste("Strain",Name,Rank,sep="_")
  strain_hist = ggplot(dataset, aes(x=strain))
  strain_hist + geom_histogram(binwidth = 0.1) + scale_x_log10()
  ggsave(file = paste(strain_title,".pdf",sep=""))
  
}