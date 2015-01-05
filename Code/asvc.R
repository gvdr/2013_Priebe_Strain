library(phytools)
mytree <-read.newick("baskerville.nex")

ids <- as.vector(read.table("//file/UsersG$/gvd16/Home/Desktop/ids.txt", quote="\""))
mytree$tip.label <- ids[,1]

mytree <- collapse.singles(mytree)

plot(mytree, cex = 0.5)

phylosig(mytree,TraitOUT[,3],method="lambda",test=TRUE)
contMap(mytree,TraitIN[,1],outline=TRUE, fsize=0.4,sig=2,legend=0)
