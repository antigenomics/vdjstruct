#!/bin/Rscript

clnum <- 8

args <- commandArgs(TRUE)
clnum <- ifelse(!is.na(args[1]), as.numeric(args[1]), clnum) 
	
dpart <- function(p1,p2)
{
	tt <- table(p1,p2)
	d <- 0
	for (i in 1:nrow(tt))
		for (j in 1:ncol(tt))
		{
			if (i<nrow(tt)) d <- d+sum(tt[i,j]*tt[(i+1):nrow(tt),j])
			if (j<ncol(tt)) d <- d+sum(tt[i,j]*tt[i,(j+1):ncol(tt)])
		}
	return(d/choose(length(p1),2))
}

pdf("generated/results/clustering_info.pdf")
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



library(gplots)
heatmap.2(mat1, main="Real clutering heatmap", dendrogram='none', trace='none', Rowv=FALSE, Colv=FALSE)
heatmap.2(mat2, main="Teoretical clutering heatmap", dendrogram='none', trace='none', Rowv=FALSE, Colv=FALSE)

plot(group1, main="Real clutering plot", xlab="Element index", ylab="Cluster number")
plot(group2, main="Teoretical clutering plot", xlab="Element index", ylab="Cluster number")

library(fpc)
sink("generated/results/clustering_log.txt", append=FALSE, split=FALSE)
sprintf("Corrected RAND:")
cluster.stats(dist1, group1, group2)$corrected.rand
sprintf("Statictics for real clustering:")
cluster.stats(dist1, group1)
sprintf("Statictics for teoretical clustering:")
cluster.stats(dist2, group2)
sink()

# Silhouette hist
library(cluster)
si <- silhouette(group1, dist2)
sils <- c()
for (i in 1:2000)
{
	group <- sample(group1)
	sils <- append(sils, mean(silhouette(group, dist2)[,3]));
}
plot(si)
hist(sils, breaks = 20, freq = FALSE, col = "lightblue", main="Silhouette distribution while shuffling and the real one", xlab="Silhouette value", ylab="Probability")
lines(density(sils), col="red")
abline(v = mean(si[,3]), col="red", lwd=6)



# Cluster num graphs
n <- 20
mydist <- rep(0,n)
for (i in 2:n)
{
  mydist[i] <- dpart(cutree(fit1, k=i), cutree(fit2, k=i))
}

plot(mydist, type='b', main="Correlation between compability of real and \nteoretical custering and clusters number\n using histogram clustering", xlab="Clusters number", ylab="Compability level", col="blue")

ks <- 2:20
WSS <- sapply(ks, FUN=function(k) {
  cluster.stats(dist2, cutree(fit1, k), cutree(fit2, k))$corrected.rand
  
})
plot(ks, WSS, type="b", main="Corrected RAND index for histogram clustering", xlab="Number of clusters", ylab="Index")



if (clnum <= 8)
{
	n <- clnum
	mydist <- rep(0,n)
	for (i in 2:n)
	{
	  mydist[i] <- dpart(kmeans(dist1, centers=i, nstart=10)$cluster, kmeans(dist2, centers=i, nstart=10)$cluster)
	}
	plot(mydist, type='b', main="Correlation between compability of real data and \nteoretical custering and clusters number\n using kmeans clustering", xlab="Number of clusters", ylab="Compability", col="blue")

	ks <- 2:clnum
	WSS <- sapply(ks, FUN=function(k) {
	  cluster.stats(dist2, kmeans(dist1, centers=k, nstart=10)$cluster, kmeans(dist2, centers=k, nstart=10)$cluster)$corrected.rand
	  
	})
	plot(ks, WSS, type="b", main="Corrected RAND index for kmeans clustering", xlab="Number of clusters", ylab="Index")
}

dev.off()
