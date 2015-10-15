#!/bin/Rscript

pdf("../../../generated/pdbcdr3/pairwise_distance_freq.pdf")
p <- "../../../generated/pdbcdr3/dist_mats/cdr3+pep"
temp_a <- list.files(path = p, pattern = "*0).txt")
matlist_a <- list()

for (i in 1:(length(temp_a)))
{
  bufer <- read.table(paste(p, temp_a[i], sep = '/'))
  if (nrow(bufer) > 1)
  {
    matlist_a[[i]] <- as.matrix(bufer[2:nrow(bufer),][,2:ncol(bufer)])
    matlist_a[[i]]
  }
}

lst <- numeric()
for (i in 1:(length(temp_a)-8))
{
  lst <- append(lst, as.numeric(matlist_a[[i]]))
}
listhist <- hist(lst, breaks = 100, 
     main = "Pairwise amino acids distance (alpha chain)\n25% cut", 
     xlab = "Distance (Angstroms)", 
     ylab = "Frequency",
     col = "lightblue")
freqs <- listhist$counts
abline(h = max(freqs)*0.25, col="blue", lwd=4)

# ==========================================================

temp_b <- list.files(path = p, pattern = "*1).txt")
matlist_b <- list()

for (i in 1:(length(temp_b)))
{
  bufer <- read.table(paste(p, temp_b[i], sep = '/'))
  if (nrow(bufer) > 1)
  {
    matlist_b[[i]] <- as.matrix(bufer[2:nrow(bufer),][,2:ncol(bufer)])
    matlist_b[[i]]
  }
}

lst <- numeric()
for (i in 1:(length(temp_b)-8))
{
  lst <- append(lst, as.numeric(matlist_b[[i]]))
}
listhist <- hist(lst, breaks = 100, 
                 main = "Pairwise amino acids distance (beta chain)\n25% cut", 
                 xlab = "Distance (Angstroms)", 
                 ylab = "Frequency",
                 col = "lightblue")
freqs <- listhist$counts
abline(h = max(freqs)*0.25, col="blue", lwd=4)

dev.off()