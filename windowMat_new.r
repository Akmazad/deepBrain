windowMatForAllchr <- function(datadir, binSize, outputPath){
    ## Load libraries
    library('RCurl')
    library('derfinder')
    library('GenomicRanges')
    
    ## first get the bam file paths
    bam_files <- list.files(path = datadir, pattern="*.sorted.bam$", recursive=TRUE, full.names=TRUE)
    ##  then make the bai paths since they are meant to be in the same directory here as the bam file
    bai_files <- paste(bam_files, ".bai", sep="")
    
    ## Load the data from disk -- for choromose 2 for example
    allChrs = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT')
    ## allChrs = c('1', 'X', 'Y', 'MT')
    

    ## summed coverage for each bin for all the samples (bam)
    for(chrID in 1:length(allChrs)){
        out = NULL
        fullCov = fullCoverage(files = bam_files, bam = bai_files, chrs = allChrs[chrID])
        start <- 1
        b <- binSize
        binID <- 1
        chrLength <- as.numeric(sapply(fullCov[chrID], nrow))
        message(paste("Chromosome ",allChrs[chrID]," has started","\n"), appendLF=FALSE)
        while (start < chrLength){
            end <- if(start + b - 1 <= chrLength) (start + b - 1) else chrLength 
            bin <- window(DataFrame(fullCov[chrID]), start, end)
            sb <- sapply(bin,sum)
            ##names(sb) <- paste(names(fullCov[chrID]), "_bin_", binID, sep="")
            out <- rbind(out,sb)
            message(paste("bins done: ",binID,", "), appendLF=FALSE)
            start <- end + 1
            binID <- binID + 1
          }
        message(paste("\n", "Chromosome ",allChrs[chrID]," has been completed","\n"), appendLF=FALSE)
        ## remove the fullCov data
        rm(fullCov)
        write.csv(out,file = paste(outputPath, "/outputMatForChr", "_", chrID, ".csv"))
      }   
}
## datadir = "/Volumes/Seagate/STAR_Output/"
## outputPath = "/home/"
## windowMatForAllchr(datadir=datadir, binSize=1000, outputPath=outputPath)
