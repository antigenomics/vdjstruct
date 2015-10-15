#!/bin/Rscript

thresh <- 8
const <- 0
savepath <- "../../../generated/pdbcdr3"

mkIntTable <- function(tablelist)
{
  alphabet <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 
                'G', 'H', 'I', 'L', 'K', 'M', 'F', 
                'P', 'S', 'T', 'W', 'Y', 'V')
  d <- matrix(0, ncol = length(alphabet), nrow = length(alphabet))
  colnames(d) <- alphabet
  rownames(d) <- alphabet
  for (i in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[i]])
    r = table[,1]
    r
    c = table[1,]
    c
    for (j in 2:length(r))
      for (k in 2:length(c))
      {
        if (as.numeric(table[j, k]) < thresh)
        {
          am1 = r[j]
          am2 = c[k]
          d[am1, am2] = d[am1, am2]+1
          d[am2, am1] = d[am2, am1]+1
        }
      }
  }
  mm = max(as.numeric(d))
  for (i in alphabet)
    for (j in alphabet)
    {
      d[i, j] = mm-d[i, j]+const
    }
  print(d)
  return(d)
}

heatTable <- function(table, name)
{
  library(gplots)
  heatmap.2(as.matrix(table), 
            dendrogram = "none", 
            Rowv = FALSE, 
            Colv = FALSE,
            main = name)
}

heatTableWithDendro <- function(table, name)
{
  library(gplots)
  #my_palette <- colorRampPalette(c("yellow", "red"))(n = 299)
  heatmap.2(as.matrix(table), #col=redgreen(255),
            #dendrogram = "none", 
            #Rowv = FALSE,
            #Colv = FALSE,
            main = name)
}

pdf(paste(savepath, "amino_interaction_matrices.pdf", sep = '/'))
p <- "../../../generated/pdbcdr3/dist_mats/cdr3+pep"

# =================== Reading matrices =====================
temp_a <- list.files(path = p, pattern = "*0).txt")
matlist_a <- list()
count <- 0
for (i in 1:(length(temp_a)))
{
  bufer <- read.table(paste(p, temp_a[i], sep = '/'))
  if (nrow(bufer) > 1)
    matlist_a[[i-count]] <- bufer
  else
    count <- count+1
}
# ==========================================================
temp_b <- list.files(path = p, pattern = "*1).txt")
matlist_b <- list()
count <- 0
for (i in 1:(length(temp_b)))
{
  bufer <- read.table(paste(p, temp_b[i], sep = '/'))
  if (nrow(bufer) > 1)
    matlist_b[[i-count]] <- bufer
  else
    count <- count+1
}
# ==========================================================

b_table <- mkIntTable(matlist_b)
a_table <- mkIntTable(matlist_a)

heatTable(a_table, "Interaction matrix for amino acids\n in alpha chain")
heatTable(b_table, "Interaction matrix for amino acids\n in beta chain")
heatTableWithDendro(b_table, "Interaction matrix for amino acids\n in beta chain (reordering allowed)")
heatTableWithDendro(a_table, "Interaction matrix for amino acids\n in alpha chain (reordering allowed)")

sink(paste(savepath, "amino_interaction_matrices.txt", sep = '/'))
print("Alpha table")
print(a_table)
print("Beta table")
print(b_table)
sink()

dev.off()