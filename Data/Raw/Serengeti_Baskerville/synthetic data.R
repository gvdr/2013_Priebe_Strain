Graph = SB_Graph
samples = 314
ranks = 6


n=vcount(Graph)
names = row.names(get.adjacency(SB_Graph))
DECO = svd(get.adjacency(Graph))
U = DECO$u
S = diag(DECO$d)
Ssqrt = structure(vapply(S, sqrt, numeric(1)),dim=dim(S))
V = DECO$v
elbows = getElbows(DECO$d[which(DECO$d > 0.001)], n=ranks)

Strains_B_lowd = array(0,dim=c(n,samples,ranks))
  
for(l in 1:ranks){
  Len_IN = Len_OUT = l+1
  
  TraitIN = V[,1:Len_IN]%*%Ssqrt[1:Len_IN,1:Len_IN]
  TraitOUT = U[,1:Len_OUT]%*%Ssqrt[1:Len_OUT,1:Len_OUT]
  p_link = TraitIN%*%t(TraitOUT)
  
  
  set.seed(42)
  for(k in 1:samples){
    Adj = p_link>=matrix(runif(n*n),nrow=n,ncol=n)
    rownames(Adj) = names
    colnames(Adj) = names
    # svds = svd(Adj)$d
    # dim = getElbows(svds[which(svds > 0.0001)],plot=F)
    Strains_B_lowd[,k,l] = strain_from_matrix(Adj,Len_OUT)
    
  }
}

par(mfrow=c(3,2))
for(i in 1:6){boxplot(Strains_B_lowd[,,i], use.cols=F)}