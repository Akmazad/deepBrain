### ba9_81.filepath <- "C:\\Users\\z3526914\\Downloads\\Brain_Prabhakar_H3K27Ac-20181113T040926Z-001\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_BA9_81_April2015_LR.csv"
### ba41_66.filepath <- "C:\\Users\\z3526914\\Downloads\\Brain_Prabhakar_H3K27Ac-20181113T040926Z-001\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_BA41_66_Mar2015_LR.csv"
### baVermis_62.filepath <- "C:\\Users\\z3526914\\Downloads\\Brain_Prabhakar_H3K27Ac-20181113T040926Z-001\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_Vermis_62_Mar2015_LR.csv"
### chrFile = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\chromosome_Length.csv"


accetylationDat <- function(ba9_81.filepath, ba41_66.filepath, baVermis_62.filepath, chrFile, binSize, overlapCutoff, outputPath){
  ### load the files
  ba9_81.dat <- read.csv(ba9_81.filepath,header = TRUE)
  ba41_66.dat <- read.csv(ba41_66.filepath,header = TRUE)
  baVermis.dat <- read.csv(baVermis_62.filepath,header = TRUE)
  
  ## read chromosome length file
  chrInfo = read.csv(file=chrFile, sep="\t", stringsAsFactors = FALSE)
  chrInd <- 1
  while(chrInd <=25 ){
    chrID <- toString(chrInfo$chrID[chrInd])
    # message(paste(chrID, " started"), appendLF=TRUE)
    chrLength <- chrInfo$chr_length[chrInd]
    # message(chrLength, appendLF=FALSE)
    start <- 1
    b <- binSize
    binID <- 1
    while (start < chrLength){
      end <- if(start + b - 1 <= chrLength) (start + b - 1) else chrLength
      
      start <- end + 1
      binID <- binID + 1
    }
    
    chrInd <- chrInd + 1
  }
  
}

