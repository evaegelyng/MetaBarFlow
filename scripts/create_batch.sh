#!/bin/bash

for FOLDER_NAME in `ls -d */`
do
  FQ_FILES=(`ls $FOLDER_NAME/*.fq`)
  echo -e "`basename ${FQ_FILES[0]}`\t`basename ${FQ_FILES[1]}`\tGGWACWRGWTGRACWNTNTAYCCYCC\tTANACYTCNGGRTGNCCRAARAAYCA\t200" > $FOLDER_NAME/batchfileDADA2.list 
  echo  "$FOLDER_NAME done"
done