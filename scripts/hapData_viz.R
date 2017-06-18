# Nandita Garud and Pleuni Pennings
# ngarud@stanford.edu and pennings@sfsu.edu
# Stanford University and San Francisco State University
# November 2014

# This code creates a visual of haplotypes and the major and minor alleles of SNPs comprising the haplotypes. This code is helpful for visualizing the missing data structure and diversity within an analysis window. 

#! /path/to/Rscript --vanilla --default-packages=utils
args = commandArgs(TRUE)
analysisWindow = read.table(args[1])
genomicData = read.table(args[2],sep=",",colClasses = "character")
outFile = args[3]
windowSize = as.numeric(args[4])
sampleSize = as.numeric(args[5])

# give each individual's haplotype data an ID ranging from 1 to the sampleSize
names(genomicData)<-c("SNPloc",paste("Strain", 1:sampleSize, sep=""))

# idenfity the center SNP of the analysis window as well as the 1st and last SNPs. 
centerSNP = analysisWindow[1,1]
firstSNP<-analysisWindow[1,2]
lastSNP<-analysisWindow[1,3]

# create a data frame for the haplotypes for this window 
haplotypes=data.frame(t(genomicData[1:windowSize,]))
names(haplotypes)<-paste("S_",genomicData$SNPloc[1:windowSize],sep="")
haplotypes<-haplotypes[2:length(haplotypes[,1]),]

# count the number of Ns per haplotypes
haplotypes$NumNs=0	
for (i in 1:length(haplotypes[,1])) {
	haplotypes$NumNs[i]<-length(which(haplotypes[i,1:windowSize]=="N"))
}

# Read in the haplotype clustering information for this analysisWindow
StrainOrder<-strsplit(as.character(analysisWindow[1,6]),",")[[1]]
NumStrainsCluster<-as.numeric(strsplit(as.character(analysisWindow[1,5]),",")[[1]])

# Count the number of singletons (NumStrainsCluster takes in a vector of the number of strains that are in haplotype clusters of 2  or greater. Need to calculate the difference between the sample size and the sum of the NumStrainsCluster). Add the singletons to NumStrainsCluster.
NumStrainsCluster<-c(NumStrainsCluster,rep(1,length(haplotypes[,1])-sum(NumStrainsCluster)))


# Parse StrainOrder by removing parentheses
for (i in 1:length(StrainOrder)){
	if (regexpr("[0-9]", StrainOrder[i])[1]!=-1) {
		StrainOrder[i]=substr(StrainOrder[i],regexpr("[0-9]", StrainOrder[i])[1],5)
		}
	else {
		StrainOrder =StrainOrder[-i]
		}
	}
	
StrainOrder<-as.numeric(StrainOrder)

# Identify which haplotypes are singletons 
StrainsNotInClusters =vector()


for (i in 1:sampleSize) {
	if (!(i %in% as.numeric(StrainOrder))) {
		StrainsNotInClusters=c(StrainsNotInClusters,i)
		}
	}

# combine all strains, either in clusters or not in clusters. 
StrainOrder<-c(StrainOrder,StrainsNotInClusters)

#Create a vector to store the cluster number (rank order of largest to smallest cluster) that each strain belongs to. Therefore, if strain 19 belongs to cluster 1, enter '1' in the 19th cell of the vector. 
haplotypes$Cluster<-0; l=1

for (clu in 1:length(NumStrainsCluster)){
	k=NumStrainsCluster[clu]
	for (j in l:(k+l-1)){
		strain=StrainOrder[j]
		haplotypes$Cluster[as.numeric(StrainOrder[j])]<-clu
		}
	l=j+1
	}


###Create haplotypes of zeros and ones. The most common haplotype will receive zeros for all SNPs comprising it, and any subsequent haplotype will receive a 1 instead of a 0 for any SNP with an alternative state from the most common haplotype.

# first find a common strain with the least amount of Ns
as.numeric(StrainOrder[1:NumStrainsCluster[1]])->commonstrains # these strains are in the most common haplotype
numNs<-vector() # vector to count the number of Ns
for (i in 1:length(commonstrains)){
	numNs=c(numNs,length(which(haplotypes[commonstrains[i],]=="N")))}	
	
StrainWMostCommonHaplo = commonstrains[which.min(numNs)] # this is the reference haplotype
ZeroOneHaplotypes=as.matrix(haplotypes)

for (i in 1:length(haplotypes[,1])){
	for (j in 1:length(haplotypes[1,])){
		if (haplotypes[i,j]=="N"){ZeroOneHaplotypes[i,j]=0} # indicate missing data with a 0
                       
			else {
                              if (haplotypes[i,j]==".")
                                  {ZeroOneHaplotypes[i,j]=8} # indicates heterozygote
                              else {
                                      if (haplotypes[i,j]==haplotypes[StrainWMostCommonHaplo,j])						{ZeroOneHaplotypes[i,j]=2} # indicate a match with the reference haplotype with a 2
                                      else {ZeroOneHaplotypes[i,j]=1}
                                  } # indicate a mismatch with the reference haplotype with a 1
                        } 
                              
	}
}
	
	
ClusterOrder<-unique(haplotypes$Cluster[sort(haplotypes$Cluster,index.return=TRUE, decreasing=TRUE)$ix])
NumStrainsClusterOrder=NumStrainsCluster[ClusterOrder]
StrainClusterOrder<-vector()
for (clu in ClusterOrder){StrainClusterOrder<-c(StrainClusterOrder,which(haplotypes$Cluster==clu))}
cols=col=c(8,9,rainbow(6),"blue","darkorange","purple","pink","darkgreen",5,6)




####
#MAKE PLOT WITH GENOTYPES (AND LATER LOCATIONS)
####
wl=0; wr=50;
cols=col=c(8,9,rainbow(6),"blue","darkorange","purple","pink","darkgreen",5,6)
#make a figure
cexpos=1
#open a pdf file to make a figure
widthpng=700; 
heightpng=1000; 	
pdf(outFile,width=5,height=7)
par(mar=c(1,1,3,1))
#make empty plot
plot(1:2,1:2,col="white",ylim=c(-10.3,sampleSize),xlim=c(-3,windowSize),ylab="",xlab="",yaxt="n",xaxt="n",frame.plot=FALSE)
title(main=paste(centerSNP),cex.main = 1.8)	
#start with most common haplotype at bottom
h=0
for (i in StrainClusterOrder){	
    height=c(h-0.35,h+0.35)
    for (j in 1:windowSize){
        rect(j-1,height[1],j,height[2],density=-1,col=as.numeric(ZeroOneHaplotypes[i,j]),lwd=0., border=NA)	
    }
    h=h+1} 

#indicate where clusters are 	
#put line for region
	h=0
	for (k in 1:length(NumStrainsClusterOrder)){
		i=NumStrainsClusterOrder[k]
		abline(h=h-0.5,col=3,lwd=1)
		rect(-2,h-0.5,-1,h+i-0.5,col=cols[floor(h%%15+1)])
		text(-15,h+0.5*i-0.5,ClusterOrder[k],cex=0.25)
		h=h+i
	}
	dev.off()
