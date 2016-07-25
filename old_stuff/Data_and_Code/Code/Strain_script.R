#read the edgelist
#the edge list has to be in "from","to"(,"weight") format
#and has to start with one row header
W_Links = as.matrix(weddell_links[,c(2,1)])
#make the graph
W_Graph = graph.edgelist(W_Links, directed=TRUE)

#if there may be unconnected nodes
#W_Graph_clean = delete.vertices(W_Graph, which(degree(W_Graph) < 1))

#choose the name
Name = "Weddell_Graph"
Graph = W_Graph

threshold <- 1.0
W_svdd = sort(svd(get.adjacency(W_Graph))$d, decreasing = TRUE)
W_svdframe = data.frame(c(1:length(W_svdd)),W_svdd)
W_svd0 <- subset(W_svdframe, value > threshold)
names(W_svdframe) = c("index","value")



ggplot() +
  geom_point(data=W_svd0, aes(x=index,y=value)) +
  geom_line(data=W_svd0, aes(x=index,y=value)) + 
  labs(title=paste(Name ,"Singular Values", sep=" "))


Traits = function(Graph,dim)
{
  Adj = get.adjacency(Graph)
  S = svd(Adj)
  Xhat1 = S$u[,1:dim] %*% diag(sqrt(S$d[1:dim]))
  Xhat2 = S$v[,1:dim] %*% diag(sqrt(S$d[1:dim]))
  return(data.frame(Xhat1,Xhat2))
}

#do the magik
#strains_C_s = NULL
for(i in 2:6){strain_plotting(Name,Graph,i)}
W_T_2 = Traits(Graph,2)
W_T_3 = Traits(Graph,3)
W_T_4 = Traits(Graph,4)
W_T_5 = Traits(Graph,5)
W_T_6 = Traits(Graph,6)
W_SD_s = data.frame(as.character(c(1:492)),SD_2$strain,SD_3$strain,SD_4$strain)
names(SD_s) = c("Node","Rank 2", "Rank 3", "Rank 4")

cols <- c("2"="blue","3"="red","4"="green")
ggplot() +
  geom_point(data=SD_2, size=3, fill="darkblue", colour="darkblue", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=SD_2, fill="blue", colour="blue", alpha=0.1, aes(x=degree, y=strain, colour="2"), method=lm, level = 0.95) +
  geom_point(data=SD_3, size=3, fill="darkred", colour="darkred", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=SD_3, fill="red", colour="red", alpha=0.1, aes(x=degree, y=strain, colour="3"), method=lm, level = 0.95) +
  geom_point(data=SD_4, size=3, fill="darkgreen", colour="darkgreen", aes(x=degree, y=strain, rank="4"), shape=21) +
  stat_smooth(data=SD_4, fill="green", colour="green", alpha=0.1, aes(x=degree, y=strain, colour="4"), method=lm, level = 0.95) + 
  labs(title="Strain vs. Node Degree") +
  scale_colour_manual(name="Rank",values=cols)

ggplot() +
  geom_point(data=SD_2, size=3, fill="darkblue", colour="darkblue", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=SD_2, fill="blue", colour="blue", alpha=0.1, aes(x=degree, y=strain, colour="2"), method=lm, formula=y ~poly(x,3), level = 0.95) +
  geom_point(data=SD_3, size=3, fill="darkred", colour="darkred", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=SD_3, fill="red", colour="red", alpha=0.1, aes(x=degree, y=strain, colour="3"), method=lm, formula=y ~poly(x,3), level = 0.95) +
  geom_point(data=SD_4, size=3, fill="darkgreen", colour="darkgreen", aes(x=degree, y=strain, rank="4"), shape=21) +
  stat_smooth(data=SD_4, fill="green", colour="green", alpha=0.1, aes(x=degree, y=strain, colour="4"), method=lm, formula=y ~poly(x,3), level = 0.95) + 
  labs(title="Strain vs. Node Degree") +
  scale_colour_manual(name="Rank",values=cols)

ggplot() +
  geom_point(data=SD_2, size=3, fill="darkblue", color="darkblue", aes(x=degree, y=strain, rank="2"), shape=21) +
  stat_smooth(data=SD_2, fill="blue", color="blue", alpha=0.1, aes(x=degree, y=strain), method=glm, formula= y ~ ns(x,3), level = 0.95) +
  geom_point(data=SD_3, size=3, fill="darkred", color="darkred", aes(x=degree, y=strain, rank="3"), shape=21) +
  stat_smooth(data=SD_3, fill="red", color="red", alpha=0.1, aes(x=degree, y=strain), method=glm, formula= y ~ ns(x,3), level = 0.95) +
  geom_point(data=SD_4, size=3, fill="darkgreen", color="darkgreen", aes(x=degree, y=strain, rank="4"), shape=21) +
  stat_smooth(data=SD_4, fill="green", color="green", alpha=0.1, aes(x=degree, y=strain), method=glm, formula= y ~ ns(x,3), level = 0.95)

ggparcoord(SD_s, columns=c(2:4), scale="globalminmax",groupColumn = "Node") + geom_line(alpha=0.5)+ xlab("Rank dimension") + ylab("Strain") + theme(legend.position="none")


ggplot() +
  geom_point(data=SB_s_names, aes(x=Species, y=X2))


SD_2_RS = regsubsets(strain ~ . , data = SD_2, nbest = 1)
plot(SD_2_RS, scale = "r2", main="Baskerville's Serengeti Strain linear regression, rank 2")

with(SD_4, invTranPlot(degree, strain, main="Baskerville's Serengeti rank 4 embedding strain"))
