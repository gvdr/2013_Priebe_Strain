
Ids <- read.csv("ids.csv", head=FALSE)
Interactions <- read.csv("interactions.csv", head=FALSE)

MatchLabels <- function(Ids,Interactions){
for(id in 1:dim(Ids)[1]){
	from_ <- Ids[id,1]
	to_ <- as.character(Ids[id,2])
	Interactions[which(Interactions$V1 == from_),1] <- to_
	Interactions[which(Interactions$V2 == from_),2] <- to_ 
}
return(Interactions)
}

Interactions <- MatchLabels(Ids,Interactions)


