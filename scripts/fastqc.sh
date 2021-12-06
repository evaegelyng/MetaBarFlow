# Script for quality checking all fastq files  
#
# Author: Eva Egelyng Sigsgaard, 6.12.2021
#
# The script runs the software FastQC (Andrews 2010) on the fastq files in each sequencing library folder, producing an html report, a text file with statistics, and a summary file
#
#!/bin/bash

for FOLDER_NAME in ./M*  # Replace "M" with the prefix of your own sequencing library folders
do
  mkdir -p "$FOLDER_NAME"/fastqc
  fastqc "$FOLDER_NAME"/*_1.fq "$FOLDER_NAME"/*_2.fq -o "$FOLDER_NAME"/fastqc
  echo  "$FOLDER_NAME done"
done
