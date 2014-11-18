#! /usr/bin/python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy
import os
import linecache

###############################


def peakFinder(inFN, outFile, threshold):
# This is the main definition of this script. In this definition, H12 peaks are identified where a peak is a tract of consecutive analysis windows lying above the threshold value specified by the user or the median H12 value if no value is specified. First all H12 values lying above the threshold are added to the dictionary H12Values, where the key corresponds to the line number of the analysis window with the H12 value and the value corresponds to the H12 value. The H12 values are then reverse sorted from highest to lowest. Then the top H12 values are successively identified, and any windows that are adjacent to the top H12 values are considered to be part of the same peak. These windows are added to the 'exclude' dictionary to avoid calling multiple peaks in the same genomic region.  The end coordinates of the peak are recorded, and any windows overlapping the peak on the left and right ends are added to the 'exclude' dictionary so that they are no considered for future peak calling.  

    # intialize dictionaries
    # windows to exclude
    exclude = {}

    # H12 values in data
    H12Values={}

    inFile      = open(inFN, 'r')

    # iterate through the file and put all the H12 values > threshold in one dictionary. Use the line number of te input file as the key
    lineNo=1
    for line in inFile:
        H12 = float(line.split('\t')[8])
        if H12 >= threshold:
            H12Values[lineNo] = H12
        else:
            exclude[lineNo]=0
        lineNo +=1

    # return a list of line numbers corresponding to the sorted H12 list
    lineNoSorted=sorted(H12Values,key=H12Values.get)
    lineNoSorted.reverse() # we want the highest H12 value first. 
    numberLines=lineNo-1 # total number of lines in the file

    # now start calling the peaks, starting with the highest H12 value
    for i in range(0, len(lineNoSorted)):

        # find the ends of the peak       
        if not exclude.has_key(lineNoSorted[i]):
            # find the right bound
            index=lineNoSorted[i]
            H12=linecache.getline(inFN,index).split('\t')[8]  

            # check that H12 is greater than threshold and that the window is not in the exclude list
            while H12>=threshold and not exclude.has_key(index) and index<=numberLines:
                H12=linecache.getline(inFN,index).split('\t')[8] 
                exclude[index]=0
                index +=1
        
            # get the right most coord
            rightMaxCoord=linecache.getline(inFN,index-1).split('\t')[2] 

            # remove subsequent windows which overlap with this peak on the right so that I can identify independent peaks
            index -=1
            subsequentLeftCoord = linecache.getline(inFN,index).split('\t')[1]
            while subsequentLeftCoord < rightMaxCoord and index <numberLines:
                exclude[index]=0
                index +=1
                subsequentLeftCoord =linecache.getline(inFN,index).split('\t')[1]
                
            #look for the left bound of the peak
            index=lineNoSorted[i]-1
            H12=linecache.getline(inFN,lineNoSorted[i]).split('\t')[8] 
        
            while H12 >=threshold and not exclude.has_key(index) and index >=1:
                H12=linecache.getline(inFN,index).split('\t')[8]
                exclude[index]=0
                index -=1

            # get the left most coord
            leftMaxCoord=linecache.getline(inFN,index+1).split('\t')[1]

            # remove subsequent windows which overlap with this peak on the left so that I can identify independent peaks
            index+=1
            subsequentRightCoord = linecache.getline(inFN,index).split('\t')[2]
            while subsequentRightCoord >leftMaxCoord and index>1:
                exclude[index]=0
                index -=1
                subsequentRightCoord = linecache.getline(inFN,index).split('\t')[2]

            # finally, exclude the window with the highest H12 value so that it is not counted again. 
            exclude[lineNoSorted[i]]=0

            # output the coordinates of this peak along with the summary statistics for the window with the highest H12 value in the peak
            newLine = linecache.getline(inFN,lineNoSorted[i]).strip('\n') + '\t' + leftMaxCoord + '\t' + rightMaxCoord + '\n'
            outFile.write(newLine)


######################
def calculateMedian(inFile):
# In this definition I calculate the median H12 value if the user fails to specify a threshold.

    tot=[] # store all the H12 values in this vector
    numLines=0
    for line in inFile:
        tot.append(float(line.split('\t')[8]))
        numLines +=1

    sortedTot=sorted(tot) # sort all the H12 values
    median=sortedTot[numLines/2]
    
    return median
            
#######################
def mkOptionParser():
    """ Defines options and returns parser """
    
    usage = """%prog  <input> 
    %prog finds H12 peaks. H12 peaks are defined as a consecutive tract of analysis windows with H12 values lying above some threshold H12 value"""

    parser = OptionParser(usage)

    parser.add_option("-t", "--threshold",      type="float",       default=-1.0,    help="Define a threshold H12 cutoff value to identify H12 peaks. A peak is defined as a consecutive tract of windows with H12 values above the threshold (default=mean H12 value in thedata)")
    parser.add_option("-o", "--outFile",        type="string",       default="-",    help="Write output to OUTFILE (default=stdout)")

    return parser



#####################


def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 1:
        parser.error("Incorrect number of arguments")


    inFN         = args[0]
    outFN         = options.outFile
    
    if inFN == '-':
        inFile = sys.stdin
    else:
        inFile      = open(inFN, 'r')
   
    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')

    threshold        = float(options.threshold)

    
    
    if threshold < 0.0: # this means that either the threshold was not specified or the user inputted a negative number. 
        threshold = calculateMedian(inFile) # the median will be used as the threshold

    
    peakFinder(inFN, outFile, threshold)


    

#run main
if __name__ == '__main__':
    main()

