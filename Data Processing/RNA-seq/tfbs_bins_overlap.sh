#!/bin/bash
bedDir=$1
fileDir=$2
#cd $fileDir
bins=$3/$4
overlapCutof=$(echo $5| bc)
#fFiles=find . -maxdepth 1 -type f
#for tfFile in $tfFiles
for tfFile in "$fileDir"/*.narrowPeak
do
  $bedDir/intersectBed -u -f $overlapCutof -a $bins.bed -b $tfFile > $tfFile.overlaps.bed
  filename="${tfFile##*/}"
  echo "$filename: [DONE]"
done

