#!/bin/Rscript

thr <- 7

pdf("generated/pdbcdr3/20_first_tables_histed(no_two_peaks)_bad.pdf")
p <- "generated/pdbcdr3/dist_mats/cdr3+pep"
temp <- list.files(path = p, pattern = "*.txt")
matlist <- list()
for (i in 1:20)
{
	bufer <- read.table(paste(p, temp[i], sep = '/'))
	matlist[[i]] <- as.numeric(as.matrix(bufer[2:nrow(bufer),][,2:ncol(bufer)]))
	hist(matlist[[i]], plot = T)
}
dev.off()
