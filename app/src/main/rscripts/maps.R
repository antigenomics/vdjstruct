
mat1 <- as.matrix(read.table("generated/results/beta_graph_real.txt", header=F))
mat2 <- as.matrix(read.table("generated/results/beta_graph_heat.txt", header=F))
jpeg("generated/results/real_heatmap.jpg")
heatmap(mat1)
dev.off()
jpeg("generated/results/clustered_heatmap.jpg")
heatmap(mat2)
dev.off()
