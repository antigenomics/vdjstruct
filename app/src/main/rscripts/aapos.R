#!/bin/Rscript


thresh <- 8
savepath <- "../../../generated/pdbcdr3"
alphabet <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 
              'G', 'H', 'I', 'L', 'K', 'M', 'F', 
              'P', 'S', 'T', 'W', 'Y', 'V')
my_palette <- colorRampPalette(c("white", "darkred"))(n = 299)
p <- "../../../generated/pdbcdr3/dist_mats/cdr3+pep"

mkCDR3PosList <- function(tablelist)
{
  aapos <- c()
  for (i in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[i]])
    pept = table[,1]
    cdr3 = table[1,]
    for (i in 2:length(cdr3))
      for (j in 2:length(pept))
      {
        if (as.numeric(table[j, i]) < thresh)
        {
          aapos <- append(aapos, i-round((length(cdr3)-1)/2+0.1))
          break;
        }
      }
  }
  print(aapos)
  return(aapos)
}

mkPeptidePosList <- function(tablelist)
{
  aapos <- c()
  for (i in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[i]])
    pept = table[,1]
    cdr3 = table[1,]
    for (j in 2:length(pept))
      for (i in 2:length(cdr3))
      {
        if (as.numeric(table[j, i]) < thresh)
        {
          aapos <- append(aapos, i-round((length(cdr3)-1)/2+0.1))
          break;
        }
      }
  }
  print(aapos)
  return(aapos)
}

mkDistTable <- function(tablelist)
{
  d <- matrix(0, ncol = 21, nrow = 21)
  colnames(d) <- -10:10
  rownames(d) <- -10:10
  for (i in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[i]])
    r = table[,1]
    c = table[1,]
    for (j in 2:length(r))
      for (k in 2:length(c))
      {
        if (as.numeric(table[j, k]) < thresh)
        {
          a = j-round((length(r)-1)/2+0.1)
          b = k-round((length(c)-1)/2+0.1)
          d[a+11, b+11] = d[a+11, b+11]+1
        }
      }
  }
  print(d)
  return(d)
}

mkHist <- function(list, name)
{
  res <- hist(list, breaks = 30, 
       main = name, 
       xlab = "Position", 
       ylab = "Frequency",
       col = "lightblue")
  return(res)
}

heatTable <- function(table, name)
{
  library(gplots)
  heatmap.2(as.matrix(table), col=my_palette,
            dendrogram = "none", 
            Rowv = FALSE, 
            Colv = FALSE,
            main = name)
}

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

a_peptlist <- mkPeptidePosList(matlist_a)
b_peptlist <- mkPeptidePosList(matlist_b)
a_cdr3list <- mkCDR3PosList(matlist_a)
b_cdr3list <- mkCDR3PosList(matlist_b)
a_dist_table <- mkDistTable(matlist_a)
b_dist_table <- mkDistTable(matlist_b)

pdf(paste(savepath, "interacting_aa_positions.pdf", sep = '/'))
mkHist(a_peptlist, "Interacting aa position relative to the center in peptide (alpha)")
mkHist(a_cdr3list, "Interacting aa position relative to the center in cdr3 (alpha)")
mkHist(b_peptlist, "Interacting aa position relative to the center in peptide (beta)")
mkHist(b_cdr3list, "Interacting aa position relative to the center in cdr3 (beta)")
heatTable(a_dist_table, "y - relative aa peptide position\nx - x - relative aa CDR3 position (alpha)")
heatTable(b_dist_table, "y - relative aa peptide position\nx - x - relative aa CDR3 position (beta)")
dev.off()