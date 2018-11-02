windowMatForAchr <- function(bam_files, bai_files, datadir, binSize, chrID, chrLength){
    ## Load libraries
    library('RCurl')
    library('derfinder')
    library('GenomicRanges')
    
    ## Load the data from disk -- for choromose 2 for example
    fullCov = fullCoverage(files = bam_files, bam = bai_files, chrs = chrID)
    ##print(fullCov)
    ##print(fullCov$chr1)

    ## summed coverage for each bin for all the samples (bam)
    out = NULL
    start <- 1
    b <- binSize
    while (start < chrLength){
        end <- if(start + b - 1 <= chrLength) (start + b - 1) else chrLength 
        ## message(paste("end position: ", end, "\n"), appendLF=FALSE)
        bin <- window(fullCov[,1],start,end)
        sb <- sapply(bin,sum)
        out <- rbind(out,sb)
        start <- end + 1
      }
}

## chrFile = "/Volumes/Seagate/STAR_Output/chromosome_Length_X_MT.csv"
## chrFile = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\chromosome_Length_X_MT.csv"
## datadir = "/Volumes/Seagate/STAR_Output/"
## datadir = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\tryData2"


windowMat <- function(datadir, chrFile, binSize, outputPath){
    ## for all the chromosome initialize the window Matrix
    mat <- NULL

    ## first get the bam file paths
    bam_files <- list.files(path = datadir, pattern="*.sorted.bam$", recursive=TRUE, full.names=TRUE)
    ##  then make the bai paths since they are meant to be in the same directory here as the bam file
    bai_files <- paste(bam_files, ".bai", sep="")

    ## read chromosome length file
    chrInfo = read.csv(file=chrFile, sep="\t", stringsAsFactors = FALSE)
    
    ## iterate through each chromosomes and get their coverage info
    chrInd <- 1
    while(chrInd <=25 ){
        chrID <- toString(chrInfo$chrID[chrInd])
        message(paste(chrID, " started"), appendLF=FALSE)
        chrLength <- chrInfo$chr_length[chrInd]
        tempMat <- windowMatForAchr(bam_files, bai_files, datadir, binSize, chrID, chrLength)
        mat <- rbind(mat, tempMat)
        chrInd <- chrInd + 1
        message(" DONE", appendLF=FALSE)
        }
    
    write.csv(out,file = paste(outputPath, "/outputMatForAll.csv"))
}
## to run this code
## windowMat(datadir=datadir, chrFile=chrFile, binSize=10000, outputPath=datadir)
