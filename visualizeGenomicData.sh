# Nandita Garud and Pleuni Pennings
# ngarud@stanford.edu and pennings@sfsu.edu
# Stanford University and San Francisco State University
# November 2014

dataFile=$1
analysisWindow=$2
coord=$3
window=$4
sampleSize=$5
outFile=$6

lineNo=`cat $dataFile | cut -f1 -d',' | grep -w $coord -n | cut -f1 -d':'`
halfWindow=`echo $(( window/2 ))`
upperLineNo=`echo $((lineNo+ halfWindow))`
doubleHalf=`echo $(( halfWindow*2 +1))`
cat $dataFile | head -$upperLineNo | tail -$doubleHalf > ${coord}_data_tmp.txt

lineNo2=`cat $analysisWindow | cut -f1 | grep -w $coord -n | cut -f1 -d':'`
cat $analysisWindow | head -$lineNo2 | tail -1 > ${coord}_H12_H2H1_tmp.txt

Rscript hapData_viz.R ${coord}_H12_H2H1_tmp.txt ${coord}_data_tmp.txt $outFile $window $sampleSize

rm -f ${coord}_data_tmp.txt 
rm -f ${coord}_H12_H2H1_tmp.txt
