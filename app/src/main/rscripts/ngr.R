#!/bin/Rscript

range <- 10
ran <- 3
thresh <- 8
lambda <- 10.0
savepath <- "../../../generated/pdbcdr3"
alphabet <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 
              'G', 'H', 'I', 'L', 'K', 'M', 'F', 
              'P', 'S', 'T', 'W', 'Y', 'V')
my_palette <- colorRampPalette(c('whitesmoke', "orange", "red"))(n = 299)
p <- "../../../generated/pdbcdr3/energy_mats"

mkCDR3PosList <- function(tablelist)
{
  aapos <- c()
  for (i in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[i]])
    pept = table[,1]
    cdr3 = table[1,]
    mins = apply(apply(table[-1, -1], c(1, 2), as.numeric), 2, min)
    for (i in 2:length(cdr3))
      for (j in 2:length(pept))
      {
        if ((as.numeric(table[j, i]) - mins[i-1] < 0.001) && (mins[i-1] < thresh))
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
    mins = apply(apply(table[-1, -1], c(1, 2), as.numeric), 1, min)
    for (i in 2:length(pept))
      for (j in 2:length(cdr3))
      {
        if ((as.numeric(table[i, j]) - mins[i-1]) < 0.001 && (mins[i-1] < thresh))
        {
          aapos <- append(aapos, i-round((length(pept)-1)/2+0.1))
          break;
        }
      }
  }
  print(aapos)
  return(aapos)
}

mkDistTable <- function(tablelist)
{
  d <- matrix(0, ncol = range*2+1, nrow = range*2+1)
  colnames(d) <- -range:range
  rownames(d) <- -range:range
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
          d[a+range+1, b+range+1] = d[a+range+1, b+range+1]+1
        }
      }
  }
  print(d)
  return(d)
}

mkRelDistTable <- function(tablelist)
{
  d <- matrix(0, ncol = 2*ran+1, nrow = ran*2+1)
  dist = matrix(NA, length(tablelist), 2*ran+1)
  colnames(d) <- -ran:ran
  rownames(d) <- -ran:ran
  colnames(dist) <- -ran:ran
  for (m in 1:length(tablelist))
  {
    table = as.matrix(tablelist[[m]])
    cuttable = apply(table[-1, -1], c(1, 2), as.numeric)
    if (min(cuttable) < thresh)
      #if (ncol(cuttable) > 13)
    {
      sum = 0 # flipping flag
      glmin = min(cuttable)
      mininds = which(cuttable == glmin, arr.ind = T)
      rmin = mininds[1, 1]
      cmin = mininds[1, 2]
      bufer <- matrix(0, ncol = 2*ran+1, nrow = ran*2+1)
      
      #if (ncol(cuttable) < 13)
      ileft = max(1, (rmin-ran))
      iright = min(nrow(cuttable),(rmin+ran))
      jdown =max(1,cmin-ran)
      jup = min(ncol(cuttable),(cmin+ran))
      
      for (i in ileft:iright)
      {
        imod = i-(rmin-ran)+1
        j = which(cuttable[i,] == min(cuttable[i, jdown:jup]))-(cmin-ran)+1
        
        bufer[imod, j] = bufer[imod, j]+1;
        #sum = sum + sign(as.numeric(colnames(d)[imod]) * as.numeric(rownames(d)[j]))
        dist[m, imod] = cuttable[i, j+(cmin-ran)-1]/glmin
      }
      #if (sum <= 0) # if cdr3 is 180 grad rotated then flip matrix from left to right
      #  bufer = t(apply(t(bufer), 2, rev))
      
      d = d + bufer
      
      #for (i in (ran-1):1)
      #  dist[m, i] = dist[m, i]/dist[m, i+1]
      #for (i in (ran+3):(2*ran+1))
      #  dist[m, i] = dist[m, i]/dist[m, i-1]
    }
  }
  
  if (F) # eliminate na rows
  {
    counter = 0
    for (m in 1:length(tablelist))
      if (all(is.na(dist[m-counter,])))
      {
        dist = dist[-(m-counter),]
        counter = counter + 1
      }
  }
  m = apply(dist, 2, function(x) mean(x, na.rm=TRUE))
  s = apply(dist, 2, function(x) sd(x, na.rm=TRUE))
  xymat = matrix(NA, 3, 2*ran+1)
  xymat[1,] = -ran:ran
  xymat[2,] = m
  xymat[3,] = s
  rownames(xymat) = c('pos', 'mn', 'dev')
  colnames(xymat) = -ran:ran
  
  print(dist)
  print(d)
  print(xymat)
  return(list("distance" = dist, "map" = d, "graph" = xymat))
}

mkMinsDistList <- function(atablelist, btablelist)
{
  dists = c()
  for (m in 1:length(atablelist))
  {
    a_table = as.matrix(atablelist[[m]])
    b_table = as.matrix(btablelist[[m]])
    a_cuttable = apply(a_table[-1, -1], c(1, 2), as.numeric)
    b_cuttable = apply(b_table[-1, -1], c(1, 2), as.numeric)
    
    #if (min(a_cuttable) < thresh && min(b_cuttable) < thresh)
    #if (ncol(cuttable) > 13)
    #{
    a_glmin = min(a_cuttable)
    b_glmin = min(b_cuttable)
    a_mininds = which(a_cuttable == a_glmin, arr.ind = T)
    b_mininds = which(b_cuttable == b_glmin, arr.ind = T)
    a_rmin = a_mininds[1, 1]
    b_rmin = b_mininds[1, 1]
    dists <- append(dists, (as.numeric(b_rmin)-as.numeric(a_rmin)))
    # }
  }
  #print(dists)
  return(dists)
}

mkHist <- function(list, name)
{
  res <- hist(list, breaks = 50, 
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
            trace='none',
            main = name)
}

mkPlot <- function(m, name)
{
  nn = data.frame(pos = c(rep(-ran:ran, length(m))), val = c(t(m)))
  library(plyr)
  cdata <- ddply(nn, 'pos', summarise,
                 N    = sum(!is.na(val)),
                 mean = mean(val, na.rm=T),
                 sd   = sd(val, na.rm=T),
                 se   = sd / sqrt(N)
  )
  print(cdata)
  library(ggplot2)
  print(ggplot(cdata, aes(x=pos, y=mean)) +
          xlab("Peptide aa position") +
          ylab("Relative distance to CDR3 (distance/global min distance)") +
          ggtitle(name) +
          geom_errorbar(aes(ymin=mean-se*3, ymax=mean+se*3), width=.1) +
          geom_line() +
          geom_point())
}

calcSumNrgs <- function(tables)
{
  nrg.list <- c()
  for (i in 1:length(tables))
  {
    table = as.matrix(tables[[i]])
    #pept = table[,1]
    #cdr3 = table[1,]
    nrg.list = append(nrg.list, sum(apply(table[-1, -1], c(1, 2), as.numeric)))
  }
  return(nrg.list)
}

calcSumSumNrgs <- function(tables.a, tables.b)
{
  nrg.list <- c()
  for (i in 1:length(tables.a))
  {
    table.a = as.matrix(tables.a[[i]])
    table.b = as.matrix(tables.b[[i]])
    #pept = table[,1]
    #cdr3 = table[1,]
    nrg.list = append(nrg.list, sum(apply(table.a[-1, -1], c(1, 2), as.numeric))+
                        sum(apply(table.b[-1, -1], c(1, 2), as.numeric)))
  }
  return(nrg.list)
}

# =================== Reading matrices =====================
temp_a <- list.files(path = p, pattern = "*0).txt")
temp_b <- list.files(path = p, pattern = "*1).txt")
matlist_a <- list()
matlist_b <- list()
count <- 0
for (i in 1:(min(length(temp_a), length(temp_b))))
{
  abufer <- read.table(paste(p, temp_a[i], sep = '/'))
  bbufer <- read.table(paste(p, temp_b[i], sep = '/'))
  if ((nrow(abufer) > 1) && (nrow(bbufer) > 1))
  {
    matlist_a[[i-count]] <- abufer
    matlist_b[[i-count]] <- bbufer
  }
  else
    count <- count+1
}
# ==========================================================
a.nrgs = calcSumNrgs(matlist_a)
mkHist(a.nrgs, 'alpha')
b.nrgs = calcSumNrgs(matlist_b)
mkHist(b.nrgs, 'beta')
all.nrgs = calcSumSumNrgs(matlist_a, matlist_b)
mkHist(all.nrgs, 'all')
#pdf(paste(savepath, ".pdf", sep = '/'))
#dev.off()
