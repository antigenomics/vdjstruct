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
p <- "../../../generated/pdbcdr3/dist_mats/cdr3+pep"

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

# =================== Compute coefficients =================
getX <- function(table, dlist, map, num)
{
  X = matrix(0, length(alphabet), length(alphabet))
  cuttable = apply(table[-1, -1], c(1, 2), as.numeric)
  if (min(cuttable) > thresh)
    return(c(rep(NA, (1+length(alphabet))*length(alphabet)/2)))
  
  sum = 0 # flipping flag
  glmin = min(cuttable)
  mininds = which(cuttable == glmin, arr.ind = T)
  rmin = mininds[1, 1]
  cmin = mininds[1, 2]
  
  rmin = rmin + num
  
  iup = max(1, (rmin-ran))
  idown = min(nrow(cuttable),(rmin+ran))
  jleft =max(1,cmin-ran)
  jright = min(ncol(cuttable),(cmin+ran))
  
  #for (i in iup:idown)
  #{
    #j = which(cuttable[i,] == min(cuttable[i, jleft:jright]))
    #imod = i-(rmin-ran)+1
    #jmod = j-(cmin-ran)+1
    #sum = sum + sign(as.numeric(c(-ran:ran)[imod]) * as.numeric(c(-ran:ran)[jmod]))
  #}
  
  cdr3 = sapply(table[1, (jleft+1):(jright+1)], as.character)
  pept = sapply(table[(iup+1):(idown+1), 1], as.character)
  
  #if (sum <= 0) # if cdr3 is 180 grad rotated then flip cdr3 from left to right
  #  cdr3 = rev(cdr3)
  
  for (i in iup:idown)
    for (j in jleft:jright)
    {
      imod = i-(rmin-ran)+1
      jmod = j-(cmin-ran)+1
      aapept = which(alphabet == pept[imod])
      aacdr3 = which(alphabet == cdr3[jmod])
      add = map[imod,jmod]/dlist[imod]
      X[aapept, aacdr3] = X[aapept, aacdr3]+add
      X[aacdr3, aapept] = X[aacdr3, aapept]+add
    }
  #colnames(X) = alphabet
  #rownames(X) = alphabet
  #print(X)
  return(as.numeric(unlist(sapply(1:nrow(X), function(x) X[x,x:nrow(X)]))))
}

getRealX <- function(table, dlist, map)
{
  return(getX(table, dlist, map, 0))
}

getFalseX <- function(table, dlist, map)
{
  return(getX(table, dlist, map, 2))
}

getXmat <- function(tablelist_a, tablelist_b, dlist_a, dlist_b, map_a, map_b)
{
  result = matrix(NA, 0, (1+length(alphabet))*length(alphabet)/2)
  y = c()
  for (i in 1:length(tablelist_a))
  {
    result = rbind(result, getRealX(tablelist_a[[i]], dlist_a, map_a)+
                     getRealX(tablelist_b[[i]], dlist_b, map_b))
    y = append(y, 1)
  }
  for (i in 1:length(tablelist_a))
  {
    result = rbind(result, getFalseX(tablelist_a[[i]], dlist_a, map_a)+
                     getFalseX(tablelist_b[[i]], dlist_b, map_b))
    y = append(y, 0)
  }
  return(list('X' = result, 'Y' = y))
}


#getXmat(matlist_a, dlist, map)

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

a_peptlist <- mkPeptidePosList(matlist_a)
b_peptlist <- mkPeptidePosList(matlist_b)
a_cdr3list <- mkCDR3PosList(matlist_a)
b_cdr3list <- mkCDR3PosList(matlist_b)
a_dist_table <- mkDistTable(matlist_a)
b_dist_table <- mkDistTable(matlist_b)
a_rel_dist_table <- mkRelDistTable(matlist_a)
b_rel_dist_table <- mkRelDistTable(matlist_b)

mins_dist <- mkMinsDistList(matlist_a, matlist_b)

sink(paste(savepath, "relative_from_most_interacting_cdr_aa_positition_matrices.txt", sep = '/'))
print("Alpha table")
print(a_rel_dist_table$map)
print("Beta table")
print(b_rel_dist_table$map)
sink()

a_rel_dist = a_rel_dist_table$distance
b_rel_dist = b_rel_dist_table$distance

pdf(paste(savepath, "interacting_aa_positions.pdf", sep = '/'))
mkHist(a_peptlist, "Interacting aa position relative to the center in peptide (alpha)")
mkHist(a_cdr3list, "Interacting aa position relative to the center in cdr3 (alpha)")
mkHist(b_peptlist, "Interacting aa position relative to the center in peptide (beta)")
mkHist(b_cdr3list, "Interacting aa position relative to the center in cdr3 (beta)")
mkHist(mins_dist, "Distance between most interacting \npeptide aa with alpha and beta")
heatTable(a_dist_table, "y - relative aa peptide position\nx - relative aa CDR3 position (alpha)")
heatTable(b_dist_table, "y - relative aa peptide position\nx - relative aa CDR3 position (beta)")
heatTable(a_rel_dist_table$map, "Relative position from the \nmost interating CDR3 aa\n y - peptide aa pos; x - CDR3 (alpha)")
heatTable(b_rel_dist_table$map, "Relative position from the \nmost interating CDR3 aa\n y - peptide aa pos; x - CDR3 (beta)")
heatTable(a_rel_dist, "Color - Distance from peptide aa in \nposition 'x' form the center to\n CDR3; x - Peptide aa position (alpha)")
heatTable(b_rel_dist, "Color - Distance from peptide aa in \nposition 'x' form the center to\n CDR3; x - Peptide aa position (beta)")
mkPlot(a_rel_dist, "Distance from peptide aa to CDR3 (alpha)")
mkPlot(b_rel_dist, "Distance from peptide aa to CDR3 (beta)")

dev.off()

library(plyr)
nn_a = data.frame(pos = c(rep(-ran:ran, length(a_rel_dist))), val = c(t(a_rel_dist)))
cdata_a <- ddply(nn_a, 'pos', summarise,
               N    = sum(!is.na(val)),
               mean = mean(val, na.rm=T),
               sd   = sd(val, na.rm=T),
               se   = sd / sqrt(N)
)
dlist_a = cdata_a$mean
map_a = a_rel_dist_table$map

nn_b = data.frame(pos = c(rep(-ran:ran, length(b_rel_dist))), val = c(t(b_rel_dist)))
cdata_b <- ddply(nn_b, 'pos', summarise,
                 N    = sum(!is.na(val)),
                 mean = mean(val, na.rm=T),
                 sd   = sd(val, na.rm=T),
                 se   = sd / sqrt(N)
)
dlist_b = cdata_b$mean
map_b = b_rel_dist_table$map

p = sample(1:1000, 210, replace=T)/10000000.

mat = getXmat(matlist_a, matlist_b, dlist_a, dlist_b, map_a, map_b)
X = mat$X
Y = mat$Y
X = na.omit(X)
m = nrow(X)

counter = 0

costfun <- function(params) {
  J = 0
  #print(counter)
  #counter = counter+1
  for (i in 1:m)
  {
    theta = t(params)
    ix = X[i,]
    mult = as.numeric(theta%*%ix)
    add = Y[i]*log(mult) + (1 - Y[i])*log(1 - mult)
    J = J - add
  }
  J = J + lambda/2.*as.numeric(theta%*%t(theta))
  J = J/m
  return(J)
}

grad <- function(params) { ## Gradient of 'fr'
  buf = c()
  for (i in 1:m)
  {
    theta = t(params)
    ix = X[i,]
    mult = as.numeric(theta%*%ix)
    buf = append(buf, 1./(1. + exp(-mult)))
  }
  return(t(X)%*%buf + lambda*params)
}

#print(costfun(p))
#grad(p)

opt = optim(p, costfun, grad, method = "CG")


