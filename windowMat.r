## datadir = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\tryData2"


windowMat <- function(datadir, binSize, chromosome){

## GOAL: TRY TO RUN WITHOUT ERROR

## Load libraries
library('RCurl')
library('derfinder')
library('GenomicRanges')

## Determine the BAM files to use (from local file)

## get the bam file paths
bam_files <- list.files(path="/Volumes/Seagate/STAR_Output/", pattern="*.bam$|*.bai$", recursive=TRUE, full.names=TRUE)
##  get the make the bai paths
bai_files <- paste(files,".bai",sep="")

## Load the data from disk -- for choromose 2 for example
fullCov = fullCoverage(files = bam_files, bam=bai_files, chrs = chromosome)


## summed coverage for each bin for all the samples (bam)
out=NULL
i <- 1
b <- binSize
while (i < 243190){
    start <- i
    end <- i + b
    bin <- window(fullCov$chr2,start,end)
    sb <- sapply(bin,sum)
    out <- rbind(out,sb)
    i <- end + 1
  }
## write.csv(out,file = "outputMat.csv")
}
