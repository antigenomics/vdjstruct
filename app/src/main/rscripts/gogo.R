#!/bin/Rscript

thresh <- 8
const <- 0
savepath <- "../../../generated/pdbcdr3"
alphabet <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 
              'G', 'H', 'I', 'L', 'K', 'M', 'F', 
              'P', 'S', 'T', 'W', 'Y', 'V')
my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 299)

mkIntTable <- function(tablelist)
{
  d <- matrix(0, ncol = length(alphabet), nrow = length(alphabet))
  colnames(d) <- alphabet
  rownames(d) <- alphabet
  pept <- c(rep(0,length(alphabet)))
  cdr3 <- c(rep(0,length(alphabet)))
  names(pept) <- alphabet
  names(cdr3) <- alphabet
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
          #d[am2, am1] = d[am2, am1]+1
        }
      }
  }
  #mm = max(as.numeric(d))
  #for (i in alphabet)
  #  for (j in alphabet)
  #  {
  #    d[i, j] = mm-d[i, j]+const
  #  }
  print(d)
  return(d)
}

mkNormTable <- function(table, freqtable)
{
  for (i in 1:length(alphabet))
    for (j in 1:length(alphabet))
    {
      p = freqtable$peptall[i]
      c = freqtable$cdr3all[j]
      if(p != 0 && c != 0)
        table[alphabet[i], alphabet[j]] = 
          table[alphabet[i], alphabet[j]]/p/c*10000
       # print(freqtable$pept[i])
    }
  return(table)
}

mkFreqTable <- function(tablelist)
{
  pept <- c(rep(0,length(alphabet)))
  cdr3 <- c(rep(0,length(alphabet)))
  peptall <- c(rep(0,length(alphabet)))
  cdr3all <- c(rep(0,length(alphabet)))
  names(pept) <- alphabet
  names(cdr3) <- alphabet
  names(peptall) <- alphabet
  names(cdr3all) <- alphabet
  for (i in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[i]])
    r = table[,1]
    r
    c = table[1,]
    c
    for (j in 2:length(r))
    {
      for (k in 2:length(c))
        if (as.numeric(table[j, k]) < thresh)
        {
          pept[r[j]] = pept[r[j]] + 1
          break
        }
      peptall[r[j]] = peptall[r[j]] + 1
    }
    for (k in 2:length(c))
    {
      for (j in 2:length(r))
        if (as.numeric(table[j, k]) < thresh)
        {
          cdr3[c[k]] = cdr3[c[k]] + 1
          break
        }
      cdr3all[c[k]] = cdr3all[c[k]] + 1
    }
  }
  d <- data.frame(pept, cdr3, peptall, cdr3all)
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

heatTableWithDendro <- function(table, name)
{
  library(gplots)
  heatmap.2(as.matrix(table), col=my_palette,
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
b_freqtable <- mkFreqTable(matlist_b)
a_freqtable <- mkFreqTable(matlist_a)
#print(names(b_freqtable$pept))
#for (i in 1:length(alphabet))
#  print(b_freqtable$pept[i])
a <- mkNormTable(a_table, a_freqtable)
b <- mkNormTable(b_table, a_freqtable)
print(a)
print(b)

heatTable(a_table, "Interaction matrix for amino acids\n in alpha chain \n(peptide - vertical)")
heatTable(b_table, "Interaction matrix for amino acids\n in beta chain \n(peptide - vertical)")
heatTableWithDendro(b_table, "Interaction matrix for amino acids\n in beta chain (reordering allowed)\n(peptide - vertical)")
heatTableWithDendro(a_table, "Interaction matrix for amino acids\n in alpha chain (reordering allowed)\n(peptide - vertical)")

heatTable(a, "Normed matrix in alpha chain \n(peptide - vertical)")
heatTable(b, "Normed matrix in beta chain \n(peptide - vertical)")
heatTableWithDendro(b, "Normed matrix in beta chain\n (reordering allowed)\n (peptide - vertical)")
heatTableWithDendro(a, "Normed matrix in alpha chain\n (reordering allowed)\n (peptide - vertical)")

sink(paste(savepath, "amino_interaction_matrices.txt", sep = '/'))
print("Alpha table")
print(a_table)
print("Beta table")
print(b_table)
sink()

sink(paste(savepath, "aa_frequency_in_interacting_pairs.txt", sep = '/'))
print("The number shows frequency with which aa appears in interacting pairs in peptide and cdr3 respectively and overall quantity of each aa (not only interacting)")
print("Alpha table")
print(a_freqtable)
print("Beta table")
print(b_freqtable)
sink()

dev.off()