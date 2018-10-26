# deepC
sample text
```R
## datadir = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\tryData2"


windowMat <- function(datadir,binSize){

## GOAL: TRY TO RUN WITHOUT ERROR

## Load libraries
library('derfinder')
library('GenomicRanges')

## Determine the BAM files to use (from local file)
files <- rawFiles(datadir = datadir,samplepatt = "*.bam$",fileterm = NULL)

## Load the data from disk -- for choromose 2 for example
fullCov = fullCoverage(files = files,chrs = '2')


## Get the region matrix of Expressed Regions (ERs):: It has three types of output: 1) regions, bpCoverage, and coverageMatrix
regionMat <- regionMatrix(fullCov, cutoff = 30, L = 76, verbose = FALSE)

## regions output as GRanges object
regionMat$chr2$regions

## Base-level coverage matrices for each of the regions
## Useful for plotting
lapply(regionMat$chr2$bpCoverage[1:2], head, n = 2)

## Dimensions of the coverage matrix
dim(regionMat$chr2$coverageMatrix)

## Coverage for each region. This matrix can then be used with limma or other pkgs
head(regionMat$chr2$coverageMatrix)


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
write.csv(out,file = "outputMat.csv",sep = "\t")
}

```R
