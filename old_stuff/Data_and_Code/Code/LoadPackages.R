load_packages <- function(parallel){
  if(parallel){
    library('doParallel')
    library('foreach')
    #start a cluster
    cl <- makeCluster(detectCores()-1) #detect the number of available cores and leave one out
    registerDoParallel(cl) #register it
  }
  
  #load all the other packages we need
  c("gridExtra",
    "magrittr",
    "phytools",
    "geiger",
    "picante",
    'plyr',
    "igraph",
    "ggplot2",
    "scales",
    "viridis"
  ) -> list.packages
  invisible(lapply(list.packages,library, character.only = TRUE))
}