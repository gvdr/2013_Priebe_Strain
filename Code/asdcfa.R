Trait <- Trait[order(rownames(Trait)),]
T_d <- dist(Trait, method = "euclidean")
B <- hclust(T_d, method="ward")
plot(B)
groups <- cutree(B, k=14)
rect.hclust(B, k=14, border="red")
cont <- table(groups,traits_cluster_13[,3])
matchClasses(cont,method="rowmax")
cont_out <- table(traits_cluster_13[,3],groups_out_W_M)
matchClasses(ccont,method="rowmax")

randIndex(cont_out)

clusters.kmeans <- kmeans(Trait, 14, algorithm="MacQueen")
plot(Trait, col=clusters.kmeans$cluster)
plot(T_PCA$scores[,1],T_PCA$scores[,2], col=clusters.kmeans$cluster, pch=clusters.kmeans$cluster)

T_Data <- as.data.frame(Trait)
T_PCA <- princomp(~.,data=T_Data,cor=TRUE)
plot(T_PCA$scores[,1],T_PCA$scores[,2], col=groups, pch=groups)
plot(T_PCA$scores[,1],T_PCA$scores[,2], col=traits_cluster_13[,3], pch=traits_cluster_13[,3])
