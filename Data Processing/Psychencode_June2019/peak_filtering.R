setwd("C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Training 2\\")
# setwd("/Volumes/Data1/PROJECTS/DeepLearning/Test/")
library(matrixStats)
library(data.table)
# object name: 'peaks'
filename <- "mergedPeakHeightMatrix_HumanFC"
load(paste0(filename, ".rda"))
val_th <- 0 # can be decided later
new_peaks <- ifelse(peaks>val_th, 1, 0); rm(peaks)
new_peaks <- ifelse(!is.na(new_peaks), 1, 0) # replacing NAs with 0
new_peaks <- new_peaks[which(rowSums2(new_peaks, na.rm = T) > 1),]

## just checking if any bin has < 2 samples peaks
#sum(rowSums2(new_peaks) < 2)
#hist(rowSums2(new_peaks), breaks = 100)
# save as a BED file
binInfo <- as.data.frame(do.call(rbind,strsplit(rownames(new_peaks), "_")))
final.dat <- cbind(binInfo, rownames(new_peaks), new_peaks)
colnames(final.dat) <- c("chr","start","end", "id", colnames(new_peaks))
fwrite(final.dat, paste0(filename, "_filtered.BED"), row.names = F, quote = F, sep="\t")
