#FP and PE fast
#based on code by Iain Martyn, 2012
#this version, Arne Mooers June 2015

library(caper)

############################    Main code - load two scripts at the end first
filepath1="~/Desktop/trees_1024.txt"
filepath2="~/Desktop/out_FP_1024.txt"
filepath3="~/Desktop/out_PE_1024.txt"

#trees are stored as a .txt set of newick strings, one tree per line, line end with ";"
AllTrees=scan(filepath1, what="",sep="\n",quiet=TRUE,skip=0,comment.char="#")
AllTrees=unlist(strsplit(AllTrees,"[;]"))

##get number of tips
tpc <- unlist(strsplit(AllTrees[1], "[\\(\\),;]"))
tpc=tpc[grep(":",tpc)]
tiplabels=tpc[grep(".:",tpc)]
nn=length(tiplabels)

# pre-allocate score matrices
SM_FP=matrix(0,nn,length(AllTrees))
SM_PE=matrix(0,nn,length(AllTrees))


###loop through all the trees
for (ii in 1:length(AllTrees))
{
p=AllTrees[ii]
B=readCAIC(p)

#Actually compute the scores
IS=FP_PE(B[[1]],B[[2]],B[[3]])
FP=IS[[1]]
PE=IS[[2]]

# sort them properly and fill the two matrices

fpp=as.matrix(FP[order(rownames(FP)),])
SM_FP[,ii]=fpp

pe=as.matrix(PE[order(rownames(PE)),])
SM_PE[,ii]=pe


}
###end of loop

### write the data
names=rownames(fpp)
rownames(SM_FP)=names
write.table(SM_FP,filepath2,col.names=FALSE,row.names=TRUE,quote=FALSE)

names=rownames(pe)
rownames(SM_PE)=names
write.table(SM_PE,filepath3,col.names=FALSE,row.names=TRUE,quote=FALSE)

############################    End of main code



####two necessary scripts follow

###### 1. FP & PE script combined for even more speed

FP_PE<- function (clade.matrix,edge.length,tip.label) {
# FP part
cset=rowSums(clade.matrix)*clade.matrix
lambda=edge.length*clade.matrix
tmp=lambda/cset
rm(lambda)
rm(cset)
tmp[is.na(tmp)]=0
FPP=as.matrix(colSums(tmp))
rownames(FPP)=tip.label

#PE part
PE=as.matrix(edge.length[1:dim(clade.matrix)[2]])
rownames(PE)=tip.label

IS=vector("list",2) #store two sets of isolation scores
IS[[1]]=FPP
IS[[2]]=PE

IS
}

####### 2. read CAIC fast scipt - stores a tree in useful format
readCAIC<-function(file) {
tree <- file
tpc <- unlist(strsplit(tree, "[\\(\\),;]"))
tpc=tpc[grep(":",tpc)]

# find the tip and edge labels
tiplabels=tpc[grep(".:",tpc)]
edgelabels=tpc[-grep(".:",tpc)]
edgelabels=c(edgelabels,":0")

# locate the clusters and edges
tree <- unlist(strsplit(tree, NULL))
x=which(tree=="(")
y=which(tree==")")
v=which(tree==":")
#these are the locations of the tips
w=setdiff(v,y+1) 

##Pass through string from left to right locating the paired parenthesis, thus
## allowing for easy assigment of edge to subtending tips
#initialize objects (M is the actual clade matrix, while E is vector with associated edge weights)
j=2 
k=length(w)+1
M=matrix(0,length(w)*2-1,length(w))
E=as.vector(matrix(0,1,length(w)*2-1))
x=c(x,y[length(y)]+1)
# Main Pass
while (length(x)>1)
{
if (x[j+1]<y[1])
{j=j+1
} else {
M[k,which(x[j]<w & w<y[1])]=1
E[k]=strsplit(edgelabels[k-length(w)],"[:]")[[1]][2]
k=k+1
y=y[-1]
x=x[-j]
j=j-1}
}

# Assign branch lengths and finished tip names to the tips
for (i in 1:length(w))
{M[i,i]=1
tmp=strsplit(tiplabels[i],"[:]")[[1]]
E[i]=tmp[2]
tiplabels[i]=tmp[1]}
M=list(M,as.numeric(E),tiplabels)
}

###end of scripts


