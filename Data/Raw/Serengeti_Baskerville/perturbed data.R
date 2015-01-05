steps = 2
#and has to start with one row header
SD_Links = as.matrix(read.csv("Serengeti_deVisser_Links.csv"))[,c(2,1)]
#make the graph
SD_Graph = graph.edgelist(SD_Links, directed=TRUE)

SD_adj = get.adjacency(SD_Graph)
N=dim(SD_adj)[1]

ps = seq(0, 1, by = 1/(steps-1))

Matrix_p_s = matrix(data=c(1:N),nrow=N,ncol=1)

set.seed(42)
for(k in 1:100){
p=0.01
Matrix_p = matrix(data=NA,nrow=N,ncol=N)
k=1
chances = runif(N*N)
  for(i in 1:N){
    for(j in 1:N){
      if(chances[k] < p){Matrix_p[i,j] = (SD_adj[i,j]+1)%%2}
      else{Matrix_p[i,j] = SD_adj[i,j]}
      k=k+1
    }
  }

    Matrix_p_s = cbind(Matrix_p_s,strain_from_matrix(Matrix_p,3))
  
}