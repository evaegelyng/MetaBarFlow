# Author: This script is based on the script C_Processing_with_DADA2.Rmd by Tobias G. Frøslev (see Frøslev et al. 2017).

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
