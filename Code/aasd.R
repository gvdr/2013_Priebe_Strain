#load the library in this order: there is a masking
library('tcltk') #for the plotting
library('shapes') #for the procrustes' bed
library('igraph') #for the graph stuff (load it after shapes as it masks V)
library('ppls')
library('ggplot2')
library('car')
library('leaps')

strain_from_matrix <- function(Matrix,Rank){
  #########################  ASGE STUFF   #########################
  Size = dim(Matrix)[1]
  #define the function computing the parameter matrix
  ASFWE = function(A,dim)
  {
    S = svd(A)
    Xhat1 = S$u[,1:dim] %*% diag(sqrt(S$d[1:dim]))
    Xhat2 = S$v[,1:dim] %*% diag(sqrt(S$d[1:dim]))
    Xhat = cbind(Xhat1, Xhat2)
  }
  
  #compute the rank Rank paramater matrix for the original fw 
  Xhat = ASFWE(Matrix,Rank)
  #now calculate the strain (via brute force)
  Vars = seq(1,2*Rank)
  strain = NULL
  for(i in 1:Size)
  {
    Xhat_i = ASFWE(Matrix[-i,-i],Rank)
    strain[i] = procOPA( Xhat_i[,Vars] , Xhat[-i,Vars] ,scale=T,reflect=T)$OSS
  }
  
  #########################  DATASET  #########################
  
  #write a csv with all the data
  
  return(strain)
}

strain_data <- function(Graph,Rank){
  
  ######################### GRAPH STUFF  #########################
  
  #basic network analysis
  Size = length(V(Graph))
  #network nodes statistics
  out_degree <- degree(Graph, mode="out")
  in_degree <- degree(Graph, mode="in")
  degree <- degree(Graph, mode="all")
  
  
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
  
  #########################  DATASET  #########################
  
  #write a csv with all the data
  dataset <- data.frame(strain,out_degree,in_degree,degree)
  
  return(dataset)
}