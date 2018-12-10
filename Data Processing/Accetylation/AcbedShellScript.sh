#!/bin/bash
#first argument which is the Bedtools bin directory
bedDir=$1	
#second argument which is the working directory
cd $2	
#third argument which is the chromosomal bed file name
bins=$3
#fourth argument which is overlap cut-off
overlapCutof=$(echo $4| bc)
#echo $overlapCutof
#fifth, sixth, and seventh arguments are b19, b41, and baVermis bedfiles
for features in $5 $6 $7
do
	$bedDir/intersectBed -u -f $overlapCutof -a $bins.bed -b $features.bed > $features.overlaps.bed
done