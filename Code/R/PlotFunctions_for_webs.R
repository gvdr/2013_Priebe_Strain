require("shapes") # It conflicts with igraph, you want it in this position
require("igraph") # To handle the webs objects and plot them
require("magrittr") # Provides pipes ( %>% )
require("NetIndices") # To compute ecological indices
require("dplyr") # Filter, select, mutate, ... data sets
 
require("reshape2") # To organise and 'shape' our data sets
require("ggplot2") # Nice plots
require("viridis")
options(stringsAsFactors = FALSE) # Ask me why



compute_strain <- function(Matrix,Rank, mode = c("total", "inward", "outward")){
  Matrix %>% dim %>% .[1] -> Size
  if(mode %>% equals("total")){
    Vars <- seq(1,2*Rank)
    asfwe <- ASFWE
  } else if(mode %>% equals("inward")){
    Vars <- seq(1,Rank)
    asfwe <- ASFWE_in
  } else if(mode %>% equals("outward")){
    Vars <- seq(1,Rank)
    asfwe <- ASFWE_out
  } else {
    stop("mode must be one of total, inward or outward")
  }
  Matrix %>% asfwe(Rank) -> Xhat
  Size %>% numeric -> strain
  for(i in 1:Size){
    Matrix %>% .[-i,-i] %>%
      asfwe(Rank) %>%
      .[,Vars] %>%
{procOPA(.,Xhat[-i,Vars], scale=T, reflect=T)$OSS} ->
  strain[i]
  }
return(strain)
}


ASFWE <- function(A,dim) {
  if(dim %>% equals(1)) {
    cbind(
      A %>% svd %$%
        "%*%"(u %>% .[,1],
              d %>% .[1] %>% sqrt) %>%
        as.matrix
      ,
      A %>% svd %$%
        "%*%"(v%>% .[,1],
              d %>% .[1] %>% sqrt) %>%
        as.matrix
    )  %>%
      return
  } else {
    cbind(
      A %>% svd %$%
        "%*%"(u %>% .[,1:dim],
              d %>% .[1:dim] %>% sqrt %>% diag)
      ,
      A %>% svd %$%
        "%*%"(v %>% .[,1:dim],
              d %>% .[1:dim] %>% sqrt %>% diag)
    )  %>%
      return
  }
}

ASFWE_out <- function(A,dim) {
  if(dim %>% equals(1)) {
    A %>% svd %$%
      "%*%"(u %>% .[,1] %>% as.matrix,
            d %>% . [1] %>% sqrt) %>%
      return
  } else { 
    A %>% svd %$%
      "%*%"(u %>% .[,1:dim],
            d %>% .[1:dim] %>% sqrt %>% diag) %>%
      return
  }
}

ASFWE_in <- function(A,dim) {
  if(dim %>% equals(1)) {
    A %>% svd %$%
      "%*%"(v %>% .[,1] %>% as.matrix,
            d %>% . [1] %>% sqrt) %>%
      return
  } else { 
    A %>% svd %$%
      "%*%"(v %>% .[,1:dim],
            d %>% .[1:dim] %>% sqrt %>% diag) %>%
      return
  }
}

get_Adj <- function(name, sep=",", header=FALSE, Links_file='interactions.csv',  Links_order=c(1,2), base='') {
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

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names)){
    return(dirs) }
  else {
    return(basename(dirs))
  }
}

Plot.foodweb <- function(Graph, dodge=0.08, size.scale=3) {
  Layer <- function(Graph){
    Graph %>% vcount -> V # How many vertices in the graph?
    Graph %>% get.adjacency(sparse = F) %>% TrophInd -> Trophs
    # We define a matrix with two columns:
    c(Trophs$OI + V %>% runif(max=dodge), # one with Omnivory Index (and a bit of noise)
      Trophs$TL - 1 # and one with the Trophic Levels
    ) %>% matrix(ncol=2) %>%
      return
  }
  Graph %>% get.adjacency(sparse = FALSE) %>% TrophInd %$% OI -> Omni
  Graph %>% Layer -> Layoutted
  par(mfrow=c(1,3))
  for(modeToPlot in c("total", "inward", "outward")){
    Graph %>%
      get.adjacency(type="both",sparse = FALSE) %>%
      as.matrix %>%
      compute_strain(3, mode = modeToPlot) %>%
      rank(ties="min") ->
      Strains
    Strains %>% max %>% viridis -> col_palette
    Strains %>% col_palette[.] -> igraph::V(Graph)$color
  Graph %>% plot.igraph(layout= Layoutted,
                        vertex.label=NA,
                        vertex.size= degree(.) %>% sqrt * size.scale,
                        edge.arrow.size=.2,
                        edge.width=.2,
                        edge.curved = TRUE,
                        edge.loop.angle = -1.5)
  }
}