```r
setwd("C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Training 2\\")
# setwd("/Volumes/Data1/PROJECTS/DeepLearning/Test/")
library(matrixStats)
library(data.table)
library(dplyr)

# object name: 'peaks'
filename <- "mergedPeakHeightMatrix_HumanFC"
load(paste0(filename, ".rda"))
val_th <- 0 # can be decided later

# get meta-data to select columns (HumanFC: control samples)
manifest <- fread("SYNAPSE_METADATA_MANIFEST.tsv", encoding = "UTF-8") %>% 
  dplyr::select(c("individualID","specimenID")) %>% 
  as.data.frame() 
meta.data <- fread("R01MH105472_HumanFC_Metadata_Clinical_agecensored.csv", encoding = "UTF-8") %>% 
  dplyr::select(c("individualID","daignosis")) %>% 
  as.data.frame()
control.specimens <- dplyr::inner_join(manifest, meta.data, by = "individualID") %>% 
  dplyr::filter(daignosis == "Control") %>% 
  dplyr::select("specimenID") %>% unlist(use.names = F)

# only apply the selected patients
new_peaks <- peaks %>% 
  as.data.frame() %>% 
  dplyr::select(control.specimens %>% 
                  as.vector()
                ) %>% 
  as.data.frame()

new_peaks <- ifelse(new_peaks>val_th, 1, 0); rm(peaks)
new_peaks <- ifelse(!is.na(new_peaks), new_peaks, 0) # replacing NAs with 0
perc_th = ncol(new_peaks) * 0.5 # at least 50% of the samples should have at least non-zero signal
peaks_1 <- which(rowSums2(new_peaks, na.rm = T) > perc_th)

dat1 = data.frame(peak_coordinate = rownames(new_peaks)[peaks_1],
                 label = "1")
dat0 = data.frame(peak_coordinate = rownames(new_peaks)[-peaks_1],
                  label = "0")

dat = rbind(dat1, dat0)
rownames(dat) = dat[,1]
dat = dat[-1]

binInfo <- as.data.frame(do.call(rbind,strsplit(rownames(dat), "_")))
final.dat <- cbind(binInfo, rownames(dat), dat)
colnames(final.dat) <- c("chr","start","end", "id", colnames(dat))
fwrite(final.dat, paste0(filename, "_filtered_signle_label.BED"), row.names = F, quote = F, sep="\t")
```
