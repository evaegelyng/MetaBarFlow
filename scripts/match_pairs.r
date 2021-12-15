# Author: This script is based on the script C_Processing_with_DADA2.Rmd by Tobias G. Frøslev (see Frøslev et al. 2017).
# The script filters reads, so that among other things, only reads with a matching paired read are retained.
# It has been modified by Eva, so both the forward and reverse files are now checked for empty files, and in case either or both files are empty, 
# the file sizes of each are printed

# Remember to check whether you would like to change any of the options of the function fastqPairedFilter. 
# See option descriptions from https://rdrr.io/bioc/dada2/man/fastqPairedFilter.html below:
# minLen: (Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after all other trimming and truncation.
# maxN: (Optional). Default 0. After truncation, sequences with more than maxN Ns will be discarded. Note that dada currently does not allow Ns.
# maxEE: (Optional). Default Inf (no EE filtering). After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
# truncQ:(Optional). Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.

args = commandArgs(trailingOnly=TRUE)

inputFiles = unlist( strsplit(args[1],",") )
outputFiles = unlist( strsplit(args[2],",") )
print("inputs")
print(inputFiles)
print("outputs")
print(outputFiles)

library(dada2) 

start.time <- Sys.time()

#Matching the "sense" reads.
if ( file.info(inputFiles[3])$size == 0 || file.info(inputFiles[4])$size == 0  ) {
  print(paste(inputFiles[3], " has size ", file.info(inputFiles[3])$size))
  print(paste(inputFiles[4], " has size ", file.info(inputFiles[4])$size)) } else {
   print(paste("processing ", inputFiles[3], "and ", inputFiles[4]))
   fastqPairedFilter(c(inputFiles[3],inputFiles[4]), c(outputFiles[3],outputFiles[4]), minLen=50, maxN=0, maxEE=2, truncQ=2,
                     matchIDs=TRUE)
   print("done")
}


#Matching the "antisense" reads.
if ( file.info(inputFiles[1])$size == 0 || file.info(inputFiles[2])$size == 0  ) {
  print(paste(inputFiles[1], " has size ", file.info(inputFiles[1])$size))
  print(paste(inputFiles[2], " has size ", file.info(inputFiles[2])$size)) } else {
   print(paste("processing ", inputFiles[1], "and ", inputFiles[2]))
   fastqPairedFilter(c(inputFiles[1],inputFiles[2]), c(outputFiles[1],outputFiles[2]), minLen=50, maxN=0, maxEE=2, truncQ=2,
                     matchIDs=TRUE)
   print("done")
}
