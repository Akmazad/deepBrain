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

# ------------- Read random Bins of given size
bin_outputFile = paste0("hg19_bins_", binSize,"bp","_rand")

# ---------------------------------------------- EpiMap
# --------------- IntersectBED with bins
EpiMap_inputFile = "mergedPeakHeightMatrix_EpiMap_filtered"
# Note: "EpiMap_inputFile" should exist within the "baseDir"
EpiMap_outputFile = paste0(EpiMap_inputFile,"_randBins.overlaps")
message("Intersect BED EpiMap:",appendLF=F)
system2('intersectBed', 
        paste('-wao -f 0.05 -a', paste0(baseDir,bin_outputFile,".bed"), '-b', paste0(baseDir,EpiMap_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,EpiMap_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)


