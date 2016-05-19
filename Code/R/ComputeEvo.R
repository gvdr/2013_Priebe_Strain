
name %>%
     paste(.,"/tree.phy",sep="") %>%
     read.newick %>%
     collapse.singles -> Tree
name %>%
     paste(.,"/ids.csv",sep="") %>%
     read.csv -> ids_table

short_names <- as.character(ids_table[,3])
long_names <- as.character(ids_table[,2])

if(all(Tree$tip.label == long_names)){
	Adj <- Adj[short_names,short_names]
	rownames(Adj) <- colnames(Adj) <- long_names
} else {
	stop("tip labels and ids names do not match")
}


EvDist <- data.frame(matrix(ncol = 3, nrow = length(long_names))
EvDist[,1:2] <- evol.distinct(Tree,type="equal.splits")
EvDist[,3] <- evol.distinct(Tree,type="fair.proportion")[,2]

