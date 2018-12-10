#!/bin/bash
bedDir=$1	#first argument which is the Bedtools bin directory
cd $2	#second argument which is the working directory
bins=$3	#third argument which is the chromosomal bed file name
overlapCutof=$(echo $4| bc)	#fourth argument which is overlap cut-off
for features in $5 $6 $7	#fifth, sixth, and seventh arguments are b19, b41, and baVermis bedfiles
do
	$bedDir/intersectBed -u -f $overlapCutof -a $bins.bed -b $features.bed > $features.overlaps.bed
done
