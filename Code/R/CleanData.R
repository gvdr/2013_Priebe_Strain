setwd("../Serengeti_Baskerville") #to be adpated to source directory
ids_file <- "ids.csv"
ids_from <- 3
ids_to <- 2


file_name <- "interactions.txt"
target_folder <- "Serengeti_Baskerville"
target_file <- paste("../../Elaborated/interactions/",target_folder,"/interactions.csv",sep="")

Interactions <- read.csv(file_name, sep="\t", head = FALSE)
Ids <- read.csv(ids_file, head = TRUE)
Interactions <- data.frame(prey = as.character(Interactions[,1]),predator = as.character(Interactions[,2]), stringsAsFactors = FALSE)

MatchLabels <- function(Ids,Interactions,from,to){
  for(id in 1:dim(Ids)[1]){
    from_ <- as.character(Ids[id,from])
    to_ <- as.character(Ids[id,to])
    if(any(Interactions$prey == from_)){
      Interactions[which(Interactions$prey == from_),1] <- to_
    }
    if(any(Interactions$predator == from_)){
      Interactions[which(Interactions$predator == from_),2] <- to_ 
    }
  }
  return(Interactions)
}

Interactions <- MatchLabels(Ids,Interactions,ids_from,ids_to)


write.table(Interactions, file = target_file, row.names = FALSE, col.names = FALSE, sep = ",")
