# Nandita Garud and Pleuni Pennings
# ngarud@stanford.edu and pennings@sfsu.edu
# Stanford University and San Francisco State University
# November 2014

#! /path/to/Rscript --vanilla --default-packages=utils

args = commandArgs(TRUE)
peaks = read.table(args[1])
outFile = args[2]
numViz = as.numeric(args[3])
sampleSize = as.numeric(args[4])

# set up the paper and divide it into the number of haplotype spectra we are visualizing
pdf(outFile,width=4.25,height=4.86, title = "Figure2",paper="special")
s=as.vector(c(1,10))
split.screen(s)
par(omi=c(0.1,0.30,0.15,0.05)) #A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.

cex2=0.6

mtext("H12:",side=1,at=-0.8,padj=10.2, cex=cex2)
mtext("H2/H1:",side=1,at=-0.8,padj=11.7, cex=cex2)

for (x in c(1:numViz)){
screen(x)

# need to adjust the margins
par(mai=c(.1,0,.05,0.03))

# parse the haplotype distribution field of the peaks output.  
a1 <- as.character(peaks[x,5]) 
a2 <- strsplit(a1, ",") 
a3 <- unlist(a2) 
hapDist <- as.numeric(a3) 

# since the H12 output shows the frequency spectrum for haplotype groups with at least two individuals, I need to count the number of haplotypes that are singletons. These singletons will be colored grey. 
NumOnes=sampleSize-sum(hapDist)
wl=0; wr=50;

cols=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#c51b7d','#7fbc41','#f1b6da','#8c510a','#c7eae5','#35978f','#dfc27d','#d6604d')

plot(1:2,1:2,col="white",ylim=c(-1.3,sampleSize),xlim=c(-1,55),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)

# color the singletons grey
for (i in 1: NumOnes){
	height=c(i-0.35,i+0.35)
	rect(wl,height[1],wr,height[2],density=-1,col="grey",lwd=0, border='grey')}

# color the haplotype groups with >= 2 members according to the color vector above. 
colNo=1
for (j in sort(hapDist)){
	height=c(height[2]+0.3,height[2]+j)
	rect(wl,height[1],wr,height[2],density=-1,col=cols[length(hapDist)-colNo+1],lwd=0, border='white')
	colNo=colNo+1
	}

# Print H12 and H2/H1 values below each bar plot. 
H12=paste(round(peaks[x,9],2))
H2H1=paste(round(peaks[x,10],2))
mtext(H12,side=1,cex=cex2,padj=-3)
mtext(H2H1,side=1,cex=cex2,padj=-1.5)
}

close.screen(all=TRUE)
dev.off() 
