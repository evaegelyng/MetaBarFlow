# Script for creating the batch files required by demultiplex.sh 
#
# Author: Eva Egelyng Sigsgaard
#
# The script produces a file (batchfileDADA2.list) in each sequencing library folder containing the fastq file names, the primer sequences, and the minimum length required for a read 
# (unaligned, i.e. forward or reverse read) after trimming of primers and tags. Replace the primer sequences and length specified in the script with those appropriate for your project. 
# If your primers contain inosine bases ("I"), these need to be replaced with "N", as the software does not recognize "I". 
#
#!/bin/bash

for FOLDER_NAME in ./M*
do
  FQ_FILES=(`ls $FOLDER_NAME/*.fq`)
  echo -e "`basename ${FQ_FILES[0]}`\t`basename ${FQ_FILES[1]}`\tGGWACWRGWTGRACWNTNTAYCCYCC\tTANACYTCNGGRTGNCCRAARAAYCA\t200" > $FOLDER_NAME/batchfileDADA2.list 
  echo  "$FOLDER_NAME done"
done
