# setwd("/Volumes/Data1/PROJECTS/Psychencode_June2019/EpiMap/ProcessedData/")
setwd("C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Training 2\\")
library(lattice)

th_select <- function(file, last_val, tfFilter=F){
  denom <- nrow(peaks)
  val_dist <- seq(1,last_val,1)
  perc_dist <- seq(0.1,1,0.05)
  mat <- matrix(data = 0,nrow = length(perc_dist), ncol = length(val_dist))
  i=0
  for(fpkm_perc_th in perc_dist){
    j=0
    for(fpkm_val_th in val_dist){
      filtered.row=which(rowMeans(peaks > fpkm_val_th) >= fpkm_perc_th)
      mat[i,j] <- length(filtered.row)/denom
      j=j+1
    }
    i=i+1
  }
  return(mat)
}

# changes
fileName <- "mergedPeakHeightMatrix_HumanFC"
dat = load(paste0(fileName,".rda"))
hist(peaks,breaks = 200)
# changes
last_val <- 20

mat <- th_select(peaks,last_val)
save(mat,file=paste0(fileName,"_peakTH.rda"))

# Plot
# pdf(paste0("valTh_VS_percTh_",file,".pdf"))
par(mar=c(3,4,2,2))
levelplot(t(mat[c(nrow(mat):1) , ]))
# dev.off() 
