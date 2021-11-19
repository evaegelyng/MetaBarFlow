#!/bin/bash

for FOLDER_NAME in `ls -d */`
do
  fastqc $FOLDER_NAME/*_1.fq $FOLDER_NAME/*_2.fq
  echo  "$FOLDER_NAME done"
done