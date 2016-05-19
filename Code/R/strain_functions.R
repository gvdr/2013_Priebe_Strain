options(stringsAsFactor = FALSE)

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
		      full.names=FALSE, ignore.case=FALSE) {
	  # use full.names=TRUE to pass to file.info
	  all <- list.files(path, pattern, all.dirs,
			               full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
    # determine whether to return full names or just dir names
    if(isTRUE(full.names))
	        return(dirs)
    else
	        return(basename(dirs))
}

get_Adj <- function(name,
		    sep=",",
		    header=FALSE,
		    Links_file='interactions.csv',
		    Links_order=c(1,2),
		    base='') {
  library('igraph') #for the graph stuff (load it after shapes as it masks V)
  Links_file <- paste(base,Links_file,sep="/")
  Links_order <- Links_order #from prey to predator
  Csv_settings <- list(sep=sep,head=header)
  #Read food web links
  Links_file %>%
    read.csv(,sep=Csv_settings$sep,header=Csv_settings$head) %>%
    as.matrix %>%
    .[,Links_order] -> Links
  ###build the graph and recover adjacency matrix###
  Graph <- graph.edgelist(Links, directed=TRUE)
  Adj <- get.adjacency(Graph)
  return(as.matrix(Adj))
}

get_Graph <- function(name,
                    sep=",",
                    header=FALSE,
                    Links_file='interactions.csv',
                    Links_order=c(1,2),
                    base='') {
  library('igraph') #for the graph stuff (load it after shapes as it masks V)
  Links_file <- paste(base,Links_file,sep="/")
  Links_order <- Links_order #from prey to predator
  Csv_settings <- list(sep=sep,head=header)
  #Read food web links
  Links_file %>%
    read.csv(,sep=Csv_settings$sep,header=Csv_settings$head) %>%
    as.matrix %>%
    .[,Links_order] -> Links
  ###build the graph and recover adjacency matrix###
  Graph <- graph.edgelist(Links, directed=TRUE)
  return(Graph)
}

ASFWE <- function(A,dim) {
  if(dim == 1) {
	S = svd(A)
	Xhat1 = as.matrix(S$u[,1]) %*% sqrt(S$d[1])
	Xhat2 = as.matrix(S$v[,1]) %*% sqrt(S$d[1])
	Xhat = cbind(Xhat1, Xhat2)
  } else { 
	S = svd(A)
	Xhat1 = S$u[,1:dim] %*% diag(sqrt(S$d[1:dim]))
	Xhat2 = S$v[,1:dim] %*% diag(sqrt(S$d[1:dim]))
	Xhat = cbind(Xhat1, Xhat2)
  }
  return(Xhat)
}

ASFWE_out <- function(A,dim) {
  if(dim == 1) {
    S = svd(A)
    Xhat1 = as.matrix(S$u[,1]) %*% sqrt(S$d[1])
    Xhat = Xhat1
  } else { 
    S = svd(A)
    Xhat1 = S$u[,1:dim] %*% diag(sqrt(S$d[1:dim]))
    Xhat = Xhat1
  }
  return(Xhat)
}

ASFWE_in <- function(A,dim) {
  if(dim == 1) {
    S = svd(A)
    Xhat2 = as.matrix(S$v[,1]) %*% sqrt(S$d[1])
    Xhat = Xhat2
  } else { 
    S = svd(A)
    Xhat2 = S$v[,1:dim] %*% diag(sqrt(S$d[1:dim]))
    Xhat = Xhat2
  }
  return(Xhat)
}

compute_mean_distances <- function(Matr, mode =c( "all", "in", "out","total", "inward", "outward"), method = "euclidean"){
	Rank <- dim(Matr)[2]/2
	if(mode %in% c("all","total")) { m_by <- Matr }
	if(mode %in% c("out","outward")) { m_by <- Matr[,1:Rank] }
	if(mode %in% c("in","inward")) { m_by <- Matr[,(Rank + 1):dim(Matr)[2]] }
	m <- as.matrix(dist(m_by, method = method))
	d <- rowMeans(m)
  return(d)
}

compute_centralities <- function(Graph, name){
  library('igraph') #for all centralities
  Size = vcount(Graph)
  dc <- degree(Graph)
  Centralities <- data.frame(Food_Web = rep(name,Size),
                             Species = igraph::V(Graph)$name,
                             DC = igraph::degree(Graph),
                             BC = igraph::betweenness(Graph),
                             CC = igraph::closeness(Graph),
                             EC = igraph::evcent(Graph, directed=TRUE)$vector,
                             IC = igraph::alpha.centrality(Graph),
                             SC = igraph::subgraph.centrality(as.undirected(Graph, mode="collapse"), diag=TRUE))
  return(Centralities)
}

compute_strain <- function(Matrix,Rank, mode = c("all", "in", "out")){
  library('shapes') #for the procrustes' bed
  Size = dim(Matrix)[1]
  if(mode == "all"){
    Vars = seq(1,2*Rank)
    asfwe <- ASFWE
  } else if(mode == "in") {
    Vars = seq(1,Rank)
    asfwe <- ASFWE_in
  } else if(mode == "out") {
    Vars = seq(1,Rank)
    asfwe <- ASFWE_out
  } else {
    stop("mode must be one of all, in or out")
  }
  Xhat = asfwe(Matrix,Rank)
  strain = NULL
  for(i in 1:Size)
  {
	Xhat_i = asfwe(Matrix[-i,-i],Rank)
		strain[i] = procOPA(Xhat_i[,Vars], Xhat[-i,Vars], scale=T, reflect=T)$OSS
  }
  return(strain)
}

get_strains_mat <- function(Mat,max_rank = 15){
  iStrain_df <- data.frame(
		Species = character(),
		Rank = numeric(),
		Strain = numeric()
		)
  N <- dim(Mat)[1]
  labels <- rownames(Mat)
  for(rank in 1:max_rank){
    rank_index <- rep(rank,N)
    istrains <- compute_strain(Mat,rank)
    iStrain_temp <- data.frame(Species = labels, Rank = rank_index, Strain = istrains)
    iStrain_df <- rbind(iStrain_df,iStrain_temp)
  }
  return(iStrain_df)
}


img_for_corr <- function(obj, col, title, lims=c(0,1)){
  image(obj, col=col, zlim= lims, main = title, asp = 1)
}

correre <- function(Mdf, value, colors, plots=TRUE, mins=TRUE){
  library("Hmisc")
  S <- dim(Mdf)[1] / 15
  x <- matrix(nrow = S, ncol = 15)
  for(i in 1:15){
    x[,i] <- Mdf[(Mdf$Rank == i),value]
  }
  Corr <- rcorr(x,type="spearman")
  if(plots == TRUE){
    img_for_corr(Corr$P, colors, paste(Mdf$Food_Web[1],Mdf$Rank[1],sep = " "))
  }
  if(mins == TRUE){
    print(paste("mean = ", mean(Corr$P, na.rm = TRUE)))
  }
  return(Corr)
}


plot_corr_df <- function(file_name, data_df, focus_value, foodws, pars){
  corrs <- vector('list', pars[1]*pars[2])
  pdf(file = file_name,width= 20, height = 10,useDingbats=F)
  par(mfrow = pars)
  i <- 1
  for(nam in foodws){
    i <- i + 1
    Temp_df <- data_df[(data_df$Food_Web == nam),c("Food_Web","Species","Rank",focus_value)]
    n_col <- dim(Temp_df)[1] / 15
    colors <- rev(dichromat(viridis(n_col), "protan"))
    corrs[[i]] <- correre(Temp_df,focus_value,colors, plots = FALSE, mins = TRUE)
  }
  dev.off()
  return(corrs)
}