windowMatForAchr <- function(bam_files, bai_files, datadir, binSize, chrID, chrLength){
    ## Load libraries
    library('RCurl')
    library('derfinder')
    library('GenomicRanges')
    
    ## Load the data from disk -- for choromose 2 for example
    fullCov = fullCoverage(files = bam_files, bam = bai_files, chrs = chrID)

    ## summed coverage for each bin for all the samples (bam)
    out=NULL
    i <- 1
    b <- binSize
    while (i < chrLength){
        start <- i
        end <- i + b
        bin <- window(fullCov$chr2,start,end)
        sb <- sapply(bin,sum)
        out <- rbind(out,sb)
        i <- end + 1
      }
}

## chrFile = "/Volumes/Seagate/STAR_Output/chromosome_Length.csv"
## datadir = "/Volumes/Seagate/STAR_Output/"

windowMat <- function(datadir, chrFile, binSize, outputPath){
    ## for all the chromosome initialize the window Matrix
    mat <- NULL

    ## first get the bam file paths
    bam_files <- list.files(path = datadir, pattern="*.sorted.bam$", recursive=TRUE, full.names=TRUE)
    ##  then make the bai paths since they are meant to be in the same directory here as the bam file
    bai_files <- paste(bam_files, ".bai", sep="")

    ## read chromosome length file
    read.csv(file=chrFile, sep="\t", stringsAsFactors = FALSE)
    
    ## iterate through each chromosomes and get their coverage info
    chrInd <- 1
    while(chrInd <=25 ){
        chrID <- chrInfo$chrID[chrInd]
        message(paste(chrID, " started"), appendLF=FALSE)
        chrLength <- chrInfo$chr_length[chrInd]
        tempMat <- windowMatForAchr(bam_files, bai_files, datadir, binSize, chrID, chrLength)
        mat <- rbind(mat, tempMat)
        chrInd <- chrInd + 1
        message(" DONE", appendLF=FALSE)
        }
    
    write.csv(out,file = paste(outputPath, "/outputMatForAll.csv"))
}
