# Nandita Garud
# ngarud@stanford.edu
# Stanford University
# November 2014

#! /path/to/Rscript --vanilla --default-packages=utils

args = commandArgs(TRUE)
H12Scan = read.table(args[1])
peaks = read.table(args[2])
outFile = args[3]
top_n = as.numeric(args[4])

pdf(outFile,width=9,height=6, title = "Figure 1",paper="special")
plot(H12Scan[,1],H12Scan[,9], pch=20, xlab='Position', ylab='H12', main='H12 scan')
points(peaks[1:top_n,1],peaks[1:top_n,9], col='red',pch=20, cex=1.5)
dev.off()