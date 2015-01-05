setwd("//file/UsersG$/gvd16/Home/Desktop/articles/carey_asge/a/Serengeti_Baskerville")
#read the edgelist
#the edge list has to be in "from","to"(,"weight") format
#and has to start with one row header
SB_Links = as.matrix(read.csv("Serengeti_Baskerville_Links.csv"))[,c(2,1)]
#make the graph
SB_Graph = graph.edgelist(SB_Links, directed=TRUE)

#choose the name
Name = "Serengeti_Baskerville"
Graph = SB_Graph


svdd = sort(svd(get.adjacency(Graph))$d, decreasing = TRUE)
svdd = data.frame(c(1:33),svdd[1:33])
names(svdd) = c("index","value")

#do the magik
#strains_C_s = NULL
Strain_2 = strain_data(Graph,2)
Strain_3 = strain_data(Graph,3)
Strain_4 = strain_data(Graph,4)
Strains = data.frame(as.character(c(1:vcount(Graph))),Strain_2$strain,Strain_3$strain,Strain_4$strain)
names(Strains) = c("Node","Rank 2", "Rank 3", "Rank 4")

samples = 314
ranks = 6


n=vcount(Graph)
names = row.names(get.adjacency(SB_Graph))
DECO = svd(get.adjacency(Graph))
U = DECO$u
S = diag(DECO$d)
Ssqrt = structure(vapply(S, sqrt, numeric(1)),dim=dim(S))
V = DECO$v
elbows = getElbows(DECO$d[which(DECO$d > 2)], n=ranks)

families = 10
samples = 4

Strains_B_random = array(0,dim=c(n,samples,families))

Len_IN = Len_OUT = 3

TraitIN = V[,1:Len_IN]%*%Ssqrt[1:Len_IN,1:Len_IN]
TraitOUT = U[,1:Len_OUT]%*%Ssqrt[1:Len_OUT,1:Len_OUT]

set.seed(42)

for(family in 1:families){
  family=1
  TraitIN_random = matrix(c(sample(TraitIN[,1]),sample(TraitIN[,2]),sample(TraitIN[,3])), ncol=3)
  TraitOUT_random = matrix(c(sample(TraitOUT[,1]),sample(TraitOUT[,2]),sample(TraitOUT[,3])), ncol=3)
  p_random_link = TraitIN_random%*%t(TraitOUT_random)
  
  for(k in 1:samples){
    Adj = p_random_link>=matrix(runif(n*n),nrow=n,ncol=n)
    rownames(Adj) = names
    colnames(Adj) = names
    Strains_B_random[,k,family] = strain_from_matrix(Adj,Len_OUT)
    
  }
}