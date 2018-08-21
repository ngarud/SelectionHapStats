# Nandita Garud
# ngarud@stanford.edu
# November 2014
# Stanford University
# This script calculates H12 and H2/H1 from population genomic data

import sys
from optparse import OptionParser
import copy
import csv
import linecache
import random
import time
import numpy
#import Bio

######################

def clusterSingleWindow(inFile, outFile, windowTot, distanceThreshold, numStrains, singleWindow):
# This definiton calculates H12, H2, and H1 in a single window centered around the coordinate specified by the user. If the coordinate is not in the file or within the bounds of a defineable window, then an error message is outputted.  

    # Count the number of lines in a file
    countLines =  open(inFile)
    numberLines =  len(countLines.readlines())
    window = int(windowTot)/2 # this is the number of SNPs on each side that I use as part of my window.
    lastSNP = numberLines - window

    center=-1
    # find the index of the center of the singleWindow
    for lineNo in range(window, lastSNP):
        coord=int(linecache.getline(inFile,lineNo).split(',')[0].strip('\n'))
        if coord == singleWindow:
            center=lineNo
    if center!=-1: # this means that the center has not be found
        flies = initialize(window, center, inFile)
        runAllDefs(flies, center, distanceThreshold, inFile, outFile, window)

    else:
        print "ERROR: Specified coordinate was not found in file or within a defineable window"

######################
def clusterHaplotypes(inFile, outFile, windowTot, jump, distanceThreshold, numStrains):
    # This definition iterates through the entire chromosome and calculates H12, H2, and H1 to each analysis window. First, we initialize the first analysis window. Since we cannot read in the entire chromosome into memory, we keep updating the 'flies' variable with data for the next analysis window by reading genomic data directly from the inFile. 
    
    # Count the number of lines in a file
    countLines =  open(inFile)
    numberLines =  len(countLines.readlines())
    window = int(windowTot)/2 # this is the number of SNPs on each side that I use as part of my window.

    lastSNP = numberLines - window
    jump = int(jump) # this is the number of SNPs I iterate by

    # Read the data into numStrains lists of Window SNPs long (200 on each side of the center)
    # center is the variable I will call to assign the middle of the sweep region I am looking at.

    # Store the haplotypes in a dictionary:
        
    center = window +1
    flies = initialize(window, center, inFile)
 
    runAllDefs(flies, center, distanceThreshold, inFile, outFile, window)

 
    ####### now iterate through the rest of the genome and fill only the part of the list that needs to be filled ########

    for center in range(window+1+jump, lastSNP, jump):
        if 2*window +1 >= jump: 

            # remove SNPs from the left
            for j in range(1,numStrains+1):
                del flies[j][0:jump]

            # Add SNPs to the right
            for i in range(0,jump):
                current_line_info = linecache.getline(inFile, center + window + i - jump + 1).split(',')
                for j in range(1,numStrains+1):
                    flies[j].append(current_line_info[j].strip())
                    
        else:
            # need to fix this -- will fail because the edge cases have not been taken care of.
            flies = initialize(window, center, inFile)

        runAllDefs(flies, center, distanceThreshold, inFile, outFile, window)

#######################
def runAllDefs(flies, center, distanceThreshold, inFile, outFile, window):
# This definition runs all other definitions to identify haplotypes, their clusters, and summary statistics. This is run for every analysis window. 

        # Count the haplotypes
        haps = countHaps(flies)

        # clump haplotypes that differ by some min threshold (including haps that differ by only an N:
        [haps_clumped, haps_clumped_count] = clusterDiffs(haps, distanceThreshold)

        # find all clusters with at least three haplotypes
        clusters = findClusters(haps_clumped)

        sizeVector = []
        keyVector = []
        if (len(clusters.keys()) == 0):
            for key in haps_clumped.iterkeys():        
                sizeVector.append(1)
                keyVector.append(key)
        else:
            [keyVector, sizeVector] = sortClusters(clusters,haps_clumped)

        centerCoord = linecache.getline(inFile,center).split(',')[0]
        edgeCoord1 = linecache.getline(inFile,center-window).split(',')[0]
        edgeCoord2 = linecache.getline(inFile,center+window).split(',')[0]
        absLengthWin = int(edgeCoord2)-int(edgeCoord1)
        

        printClusters(inFile, outFile, centerCoord, clusters, haps_clumped,  keyVector, sizeVector, absLengthWin, edgeCoord1, edgeCoord2)



#######################
def initialize(window, center, inFile):
# This definition intializes the flies dictionary. Flies takes in the strain number as the key and the haplotype as the value. The flies vector is populated with SNPs read directly from the inFile, using center as a marker as to which SNPs to read in. 

    flies = {}
    for i in range(1,numStrains+1):
        flies[i] = []

    # Add SNPs to the left and the right
    for i in range(0,window+1):
        #left
        for j in range(1,numStrains+1):
            flies[j].append(linecache.getline(inFile,center-window +i).split(',')[j].strip())
        

    for i in range(1,window+1):
        #right
        for j in range(1,numStrains+1):
            flies[j].append(linecache.getline(inFile,center +i).split(',')[j].strip())


    return flies
#######################
        
def countHaps(flies):
# In this definition, I will collapse all unique haplotypes into single instances in the haps dictionary (the key) and use a value to indicate the number of individuals with this haplotype. 

        # dictionary to store all haplotypes (to count max)
        haps = {}
        for j in range(1,numStrains+1):
            line = ''.join(flies[j])
            haps.setdefault(line,[]) # store in an array the line numbers corresponding to the haplotypes that comprise a cluster. Line numbers correspond to the processed data matrix (SNPs).   
            haps[line].append(j)
        
        return haps
        
#################

def clusterDiffs(haps, distanceThreshold):

   # In this definition I will cluster haplotypes that differ by some min threshold. If a haplotype matches another haplotype at all positions except for sites where there are Ns (missing data), then the haplotypes will be combined and the 'distance' between the two haplotypes will be considered 0. Only ATGC differnces between haplotypes will count towards the distance threshold. 

   
    distanceThreshold = int(distanceThreshold)
    haps_clumped = {} # stored all the clumped haplotypes in this hash. I will pass this into def findClusters later on. 
    haps_clumped_count = {} # I would like to record the number of different unique haplotypes that are clumped -- will htis help me later to distinguish ancestral haplotypes?

    #  I need to keep track of which key has been compared
    compared = {}

    # Now calculate the distance between unique clustering haplotypes
    for key1 in haps.iterkeys():
        if (key1 in compared) == False:
            compared[key1]=1
            haps_clumped[key1] = haps[key1]  # regardless of whether or not key1 matches anything, I need to include it in haps_clumped. Therefore I will initialize it with it's own array.
            haps_clumped_count[key1] = 1

            
            for key2 in haps.iterkeys():
                if ((haps[key2][0] in haps_clumped[key1]) == False) and ((key2 in compared) == False):
                    [distance, s1]= hamming_distance_clump(key1, key2, distanceThreshold)
                    
                    # If I replace an "N" in key1, I will replace the returned key1 in haps_clumped:
                    if distance == 0 and key1 != s1:
                        haps_clumped_count[s1] =  haps_clumped_count[key1]
                        haps_clumped[s1] = haps_clumped[key1]
                        del haps_clumped_count[key1]
                        del haps_clumped[key1]
                        key1 = s1
                    if distance <= distanceThreshold:
                        # The reason why this extra if statement is here is so that I do not confuse merging missing data with clumping haplotypes with a min distance threshold
                        # store into the haps_clumped threshold:
                        haps_clumped[key1] += haps[key2] # add the array for key2 to key1 array
                        haps_clumped_count[key1] += 1
                        compared[key2] = 1 # this means that I won't check this distance again since it has been clumped. 
    return [haps_clumped, haps_clumped_count]
                
        

##################

def findClusters(haps):
# This definition identifies haplotypes present in the sample in at least 2 individuals.

        n_min1=2
        n_min2=2
        # find the top clusters comprised of at least n_min members 
        clusters = {}
    
        # flag for first cluster > n_min1 found
        n_min1_found = False

        for key in haps.iterkeys():
            if len(haps[key]) > int(n_min1)-1:
                clusters[key] = [] # Store the top clusters in this dictionary
                n_min1_found = True

        if n_min1_found == True:
            for key in haps.iterkeys():
                if (len(haps[key]) > int(n_min2)-1 and len(haps[key]) < int(n_min1) ):

                    clusters[key] = []

        return clusters
####################
def sortClusters(clusters, haps):
# this definition sorts haplotype clusters in reverse order from largest to smallest. This sorting will help in the computation of haplotype homozygosity statistics H12, H2, and H1. 

    
        # First order the keys for each cluster from largest to smallest:
        # Put everything in vectors that I can sort

        keyVector = []
        sizeVector = []
        for key in clusters.iterkeys():
            keyVector.append(key)
            sizeVector.append(len(haps[key]))

        # now sort using bubble sort (need to sort in place):
        swapped = True
        while swapped == True:
            swapped = False
            for i in range(0, len(sizeVector)-1):
                if sizeVector[i] < sizeVector[i+1]:
                    tmpSize = sizeVector[i]
                    sizeVector[i] = sizeVector[i+1]
                    sizeVector[i+1] = tmpSize
                    tmpKey = keyVector[i]
                    keyVector[i] = keyVector[i+1]
                    keyVector[i+1]=tmpKey
                    swapped = True
        
        return [keyVector, sizeVector]


######################
def printClusters(inFile, outFile, centerCoord, clusters, haps,  keyVector, sizeVector, absLengthWin, edgeCoord1, edgeCoord2):

# this definition calculates several summary statistics from the haplotypes including K (number of unique haplotypes) and haplotype homogygosity statistics H12, H1, and H2. The outputting of all summary statistics and relevant information about the analysis window is done in this definition as well. 

        clusterSize = ''
        membersOfClusters =  ''
        hapsClumpedCount = ''

        
        for key in keyVector:
            sort_haps_vector=sorted(haps[key])
            membersOfClusters += '('
            clusterSize += str(len(haps[key])) + ','
            for k in range(0, len(sort_haps_vector)):
                membersOfClusters +=  str(sort_haps_vector[k]) + ','

            membersOfClusters += ')'
            

        # Calculate the sum of the square of the sizes of the clusters, then mulitply the first and second sizes
        H1 =0
        H2 =0
        H12 = 0
        H123 =0
        ratioH2H1 =0

        H1_vector = []
        
        for y in range(0,len(sizeVector)):
            H1 +=  (float(sizeVector[y])/float(numStrains))**2
            H1_vector.append((float(sizeVector[y])/float(numStrains))**2) 
        # Add on the singletons to H1:
        H1 += (float(numStrains)-sum(sizeVector))*(1/float(numStrains))**2
        
        
        if len(sizeVector) >0 :
            H2 = H1 - (float(sizeVector[0])/float(numStrains))**2
            
        else:
            H2 = H1 - (1/float(numStrains))**2
            
            
        if len(sizeVector) >1 :
            # calculate H12:
            H12 = ((float(sizeVector[0])+float(sizeVector[1]))/float(numStrains))**2
        elif len(sizeVector) ==1:  # If there is only 1 cluster to be found, H12 = H1  
            H12 = ((float(sizeVector[0]))/float(numStrains))**2 
        else: # if there are zero clusters:
            H12 = (2/float(numStrains))**2 -2*(1/float(numStrains))**2

        if len(sizeVector) >2:
            # finish computing H12
            for y in range(2,len(sizeVector)):
                H12 += (float(sizeVector[y])/float(numStrains))**2
            # compute H123:
            H123 = ((float(sizeVector[0])+float(sizeVector[1])+float(sizeVector[2]))/float(numStrains))**2
        
        # add on the singletons:
        H12+= (float(numStrains)-sum(sizeVector))*(1/float(numStrains))**2

        if len(sizeVector) >3:
            # finish computing H123:
            for y in range(3,len(sizeVector)):
                H123 += (float(sizeVector[y])/float(numStrains))**2

        # add on the singletons
        H123+=(float(numStrains)-sum(sizeVector))*(1/float(numStrains))**2
        
        if H1 > 0:
            ratioH2H1 = float(H2)/float(H1)


        #In the case where there are zero clusters:
        if len(clusters.keys()) == 0:
            a = str(0)
            membersOfClusters = '()'

        # print the total number of unique haplotypes found in the sample (K)
        K = len(sizeVector) + numStrains - sum(sizeVector)
        

        outFile.write(centerCoord + '\t' + edgeCoord1 + '\t' + edgeCoord2  + '\t' + str(K) +  '\t' + clusterSize + '\t' +  membersOfClusters + '\t' + str(H1) + '\t' + str(H2) + '\t'  + str(H12) + '\t' + str(ratioH2H1) + '\t' + str(H123) +  '\n')


        

#############################
def hamming_distance_clump(s1, s2, distanceThreshold):
    # this definition combines haplotypes that differ by some hamming distance (distanceThreshold as defined by the use) into a single haplotype. This is a useful feature to use if your data has a lot of sequencing errors, which may artificially inflate the number of unique haplotypes in a sample. Note that a single haplotype comprised of the union of the two haplotypes will be returned, where the ith index is replaced with the nucleotide belonging to the second haplotype being compared. 

        distance = 0
        for i in range(len(s1)):
            if s1[i] != s2[i]: 
                if (s2[i] != 'N'):
                    if (s1[i]!='N'):

                        distance += 1
                        if distance > distanceThreshold:
                            return [distance, s1]
                    else:
                        
                        s1 = s1[:i] + s2[i] + s1[i+1:] 
                        
        return [distance, s1]




######################
def mkOptionParser():
    """ Defines options and returns parser """

    usage = """%prog  <input> <number of strains> 
    %prog calculates haplotype homozygosity chromosome-wide in the window size specified. The following summary statistics are outputted: left edge coordinate of window, right edge coordinate of window, center coordinate of window, K (number of unique haplotypes), haplotype frequency spectrum (singletons are not outputted for brevity), members of haplotypeGroups, H1, H2, H12, H2/H1. """

    parser = OptionParser(usage)

    parser.add_option("-o", "--outFile",        type="string",       default="-",    help="Write output to OUTFILE (default=stdout)")
    parser.add_option("-w", "--window",        type="int",       default=400,    help="Analysis window size in terms of SNPs. Note that if the window size is an even integer, then 1 SNP will be added in order to calculate a center (default=400)")
    parser.add_option("-j", "--jump",        type="int",       default=50,    help="Distance between the centers of analysis windows in terms of SNPs (default=50)")
    parser.add_option("-d", "--distanceThreshold",      type="int",       default=0,    help="Define a threshold hamming distance such that haplotypes that are different by this amount will be grouped together as a single haplotype. This should be used in cases where the probability of sequencing errors is high (default=0, which means that any haplotype different by even 1 nucleotide will be considered its own haplotypes)")
    parser.add_option("-s", "--singleWindow",      type="int",       default=-100,    help="Instead of calculating H12 genome-wide, calculate H12 in a single window centered around the coordinate specified. The coordinate must exist in the input file in order for any output to be made")


    return parser

########################################################################

def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments")


    inFile       = args[0]
    global numStrains
    numStrains = int(args[1])

    outFN        = options.outFile
    windowTot    = options.window
    jump         = options.jump
    distanceThreshold = options.distanceThreshold
    singleWindow = options.singleWindow

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')

    if singleWindow !=-100:
        clusterSingleWindow(inFile, outFile, windowTot, distanceThreshold, numStrains, singleWindow)
    else:
        clusterHaplotypes(inFile, outFile, windowTot, jump, distanceThreshold, numStrains)


    
#run main
if __name__ == '__main__':
    main()
