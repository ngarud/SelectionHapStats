SelectionHapStats README
=================

Nandita Garud
ngarud@stanford.edu
November 2014
Petrov Lab, Stanford University

This script calculates haplotype homozygosity statistics H12 and H2/H1 from population genomic data (see Introduction below for a detailed description). Further description of H12 and H2/H1 and its application to Drosophila melanogaster population genomic data can be found in the arXived paper: Recent selective sweeps in North American Drosophila melanogaster show signatures of soft sweeps (http://arxiv.org/abs/1303.0906). 

Introduction:

Evolutionary adaptation is a process in which beneficial mutations increase in frequency in response to selective pressures. If these mutations were previously rare or absent from the population, adaptation should generate a characteristic signature in the genetic diversity around the adaptive locus, known as a selective sweep. Such selective sweeps can be distinguished into hard selective sweeps, where only a single adaptive mutation rises in frequency, or soft selective sweeps, where multiple adaptive mutations at the same locus sweep through the population simultaneously. Here we present a new statistical method that can identify both hard and soft sweeps in population genomic data. 

Our statistical test for detecting hard and soft sweeps is based on the reasoning that selective sweeps will elevate the frequencies of the genetic material surrounding the adaptive mutation, also known as a haplotype. The increase of haplotype population frequencies in both hard and soft sweeps can be captured using haplotype homozygosity. If pi is the frequency of the ith most common haplotype in a sample, and n is the number of observed haplotypes, then haplotype homozygosity is defined as H1 = Σi=1,…n pi^2. However, as we demonstrate in our paper, H1 has better ability to detect hard sweeps than soft sweeps. To have a better ability to detect hard and soft sweeps using homozygosity statistics, we developed a modified homozygosity statistic, H12 = (p1 + p2)2 + Σi>2 pi2 = H1 + 2p1p2, in which the frequencies of the first and the second most common haplotypes are combined into a single frequency.

In order to gain intuition about whether a sweep identified with H12 can be easily generated by hard sweeps versus soft sweeps under several evolutionary scenarios, we developed a new homozygosity statistic, H2/H1, where H2 = Σi>1 pi2 = H1 – p12 is haplotype homozygosity calculated using all but the most frequent haplotype. We expect H2 to be lower for hard sweeps than for soft sweeps because in a hard sweep, only one adaptive haplotype is expected to be at high frequency and the exclusion of the most common haplotype should reduce haplotype homozygosity precipitously. When the sweep is soft, however, multiple haplotypes exist at high frequency in the population and the exclusion of the most frequent haplotype should not decrease the haplotype homozygosity to the same extent. Conversely H1, the homozygosity calculated using all haplotypes, is expected to be higher for a hard sweep than for a soft sweep as we described above. The ratio H2/H1 between the two should thus increase monotonically as a sweep becomes softer, thereby offering a summary statistic that in combination with H12 can be used to test whether the observed haplotype patterns are likely to be generated by hard or soft sweeps. Note that we intend H2/H1 to be measured near the center of the sweep where H12 is the highest, otherwise further away from the sweep center mutations and recombination events will decay the haplotype signature and hard and soft sweep signatures may look indistinguishable. In order to assess if a sweep’s H12 and H2/H1 values can be easily generated under hard and soft sweeps, the H12 and H2/H1 values generated in extensive simulations under a broad range of evolutionary scenarios of sweeps of varying softness must be compared with the observed values of H12 and H2/H1.

Preprocessing steps: 
H12_H2H1.py takes in an input file consisting of all the SNPs in a population sample (see below for the format). Recommended preprocessing steps include:
1.	Remove whole strains that have very poor data quality and high amounts of missing data. Lots of missing data can result in artificially inflated H12 values. 
2.	Convert any non-ATGC site into an ‘N’ (i.e. Y should be changed to an N)
3.	Remove any closely related strains such as siblings or first-cousins 
4.	Remove any invariant sites

Running the script:
python H12_H2H1.py –help 

Required input arguments:
Input file 
Sample size

Optional input arguments (though, it is highly recommended you specify the values most appropriate for your data):
--outFile, -o: Output file name (default is stdout).
--window, -w: Analysis window size in terms of SNPs. Note that if you specify an even integer, then 1 additional SNP will be added as the center SNP (default is 400).
--jump, -j: Distance between the centers of analysis windows (in terms of SNPs) (default is 50). 
--distanceThreshold, -d: Tolerable hamming distance between unique haplotypes  (default is 0).
--singleWindow, -s: Calculate H12 in a single analysis window instead of genome-wide, centered around a specified coordinate. If the coordinate is not in the file or is too close to the ends of the files such that a definable window cannot be created, then an error is given. 

Input file format:
The input file is a line-by-line listing of all the coordinates on the chromosome where there is a polymorphism and the nucleotide states in each individual in the sample. For example:
CoordX, ind1, ind2, ind3, .. etc

A file could look like this:
Coord1,A,A,A,N,T,A,A,A,T,…etc
Coord2,A,N,G,G,G,G,G,A,G,…etc
Coord3,T,G,T,T,T,G,G,T,G,T,…etc
Etc..
Example run of the script:
python H12_H2H1.py testIn.txt 145 –o testOut.txt –w 401 –j 50 –d 0 

To calculate H12 in a single window defined at a particular coordinate:
python H12_H2H1.py testIn.txt 145 –o testOut.txt –w 401 –j 50 –d 0 –s 53034
(Here the coordinate is 53034)

Output fields:  
Field 1: Center position of the analysis window
Field 2: Left coord of analysis window 
Field 3: Right coord of analysis window
Field 4: K (number of unique haplotypes) 
Field 5: Haplotype frequency spectrum (The number of haplotypes in each haplotype group is separated by ',')
Field 6: DGRP strain numbers in each of the haplotype groups, separated by '( )' . The strain numbers correspond to the order in which the genotypes are entered in the input file (i.e. 1= ind1)
Field 7: H1 (same as regular homozygosity) 
Field 8: H2 
Field 9: H12  
Field 10: H2/H1  

Interpreting the H12 results and recommended post-processing steps:
To visualize the output, you can plot the H12 values calculated per analysis window in a program such as R (see example code below). However, there are a few recommended post-processing steps:

1.	Compare observed H12 values with a distribution of H12 values measured in simulations of realistic neutral demographic scenarios. If your observed H12 values are significantly greater than the distribution of H12 values generated under the neutral demographic model, then the neutral demographic model alone cannot explain the elevation in haplotype homozygosity in the data. 
2.	Identify extreme outliers, or “H12 peaks” in the data. To do so, run peakFinder.py, which identifies consecutive tracts of H12 analysis windows that lie above the median H12 value in the data (see README_peakFinder.docx for further details on how to run this script). 
3.	Remove regions with low recombination rates. These regions can result in artificially high H12 values. Recombination maps can either be estimated from the data, but more preferably, can be measured in an independent sample using crossing-over methods. 
4.	Check if there is any correlation between the locations of peaks with locations of inversions in the data.

Here is example R code to plot H12 values and their peaks:

Interpreting the H2/H1 results:
To gain intuition as to whether sweeps resemble hard or soft sweeps, the first thing to do is to visualize the frequency spectra of the most extreme peaks. 

Here is example R code to visualize the haplotype frequency spectra of the most extreme peaks. 

The H12 and H2/H1 statistics must be used in conjunction with one another in order to assess if a sweep’s H12 and H2/H1 values can be easily generated under hard and soft sweeps. Using approximate Bayesian computation, the H12 and H2/H1 values in data can be compared with H12 and H2/H1 values generated in simulations under a broad range of evolutionary scenarios of sweeps of varying softness. In general, both a high H12 and H2/H1 is indicative of a sweep compatible with a soft sweep, but absolute values of H12 and H2/H1 must be compared with distributions generated in simulations in order to make any contextual sense (see Garud et al. 2014 for further description on application of H2/H1 to Drosophila data).
  
