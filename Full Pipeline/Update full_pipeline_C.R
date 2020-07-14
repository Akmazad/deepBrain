# Input:  binsize, flankingLength, howManyBins (default: 1M)
#         mergedPeaks (mergedPeakHeightMatrix_HumanFC_filtered.bed, 
#         mergedPeakHeightMatrix_EpiMap_filtered.bed, 
#         Brain_CagePeaks_filtered.bed,
#         final.tf.overlaps.bed)
# Host: RNA machine/KATANA 
#       RNA:    Data dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#               Script dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#       KATANA: Data dir: /srv/scratch/z3526914/DeepBrain/Data
#               Script dir: /srv/scratch/z3526914/DeepBrain/Scripts
options(scipen=999) # prevent auto-scientific notations of numeric values

library(dplyr)
library(data.table)
source = "EpiMap"
# ------------- Set Bin-flanking configuration
binSize = 400
flanking = 2*binSize     # can be arbitrarily given
baseDir = "/srv/scratch/z3526914/DeepBrain/Data/"
dataDir = paste0(baseDir, source, "/",binSize,"_",flanking,"/")
dir.create(dataDir, recursive=T) # create the directory if doesn't exists
