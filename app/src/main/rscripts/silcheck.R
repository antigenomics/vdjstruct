
clnum <- 8

mat1 <- as.matrix(read.table("generated/results/beta_graph_real.txt", header=F))
dist1 <- as.dist(mat1)
fit1 <- hclust(dist1)
group1 <- cutree(fit1, clnum)
plot(as.dendrogram(fit1), main="Real distribution dendrogram", leaflab="none")
rect.hclust(fit1, k=clnum)

mat2 <- as.matrix(read.table("generated/results/beta_graph_heat.txt", header=F))
dist2 <- as.dist(mat2)
fit2 <- hclust(dist2)
group2 <- cutree(fit2, clnum)
plot(as.dendrogram(fit2), main="Teoretical distribution dendrogram", leaflab="none")
rect.hclust(fit2, k=clnum)

library(cluster)
si <- silhouette(group1, dist2)
sils <- c()
for (i in 1:2000)
{
	group1 <- sample(group1)
	sils <- append(sils, mean(silhouette(group1, dist2)[,3]));
}
plot(si)
hist(sils, breaks = 20, freq = FALSE, col = "lightblue")
lines(density(sils), col="red")
abline(v = mean(si[,3]), col="red", lwd=6)
