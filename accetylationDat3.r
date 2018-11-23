## preconditions: 
## 1. overlapCutoff value is between [0,1] and applied as the percentage of binSize 
## 2. window (i.e. [wStart,wEnd]) will not completely engulf the chromosomal region (i.e. [rStart,rEnd])
## 3. Hence, there will not be more than 1 overlap for each window

# ba9_81.filepath <- "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_BA9_81_April2015_LR.csv"
# ba41_66.filepath <- "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_BA41_66_Mar2015_LR.csv"
# baVermis_62.filepath <- "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_Vermis_62_Mar2015_LR.csv"
# samplefilePath = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\BrainSampleList.csv"
# chrFile = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\chromosome_Length.csv"
# outputPath = "C:\\Users\\z3526914\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data"

ba9_81.filepath <- "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_BA9_81_April2015_LR.csv"
ba41_66.filepath <- "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_BA41_66_Mar2015_LR.csv"
baVermis_62.filepath <- "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\normalized_log2_tags_Vermis_62_Mar2015_LR.csv"
samplefilePath = "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\Brain_Prabhakar_H3K27Ac\\BrainSampleList.csv"
chrFile = "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data\\chromosome_Length.csv"
outputPath = "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Data"

# ba9_81.filepath <- "/short/yr31/aa7970/azData/DeepBrain/Data/normalized_log2_tags_BA9_81_April2015_LR.csv"
# ba41_66.filepath <- "/short/yr31/aa7970/azData/DeepBrain/Data/normalized_log2_tags_BA41_66_Mar2015_LR.csv"
# baVermis_62.filepath <- "/short/yr31/aa7970/azData/DeepBrain/Data/normalized_log2_tags_Vermis_62_Mar2015_LR.csv"
# samplefilePath = "/short/yr31/aa7970/azData/DeepBrain/Data/BrainSampleList.csv"
# chrFile = "/short/yr31/aa7970/azData/DeepBrain/Data/chromosome_Length.csv"
# outputPath = "/short/yr31/aa7970/azData/DeepBrain/Data/"

accetylationDat <- function(ba9_81.filepath, ba41_66.filepath, baVermis_62.filepath, samplefilePath, chrFile, binSize, overlapCutoff, outputPath){
  ba9_81.dat <- read.csv(ba9_81.filepath,header = TRUE, stringsAsFactors = FALSE)
  ba41_66.dat <- read.csv(ba41_66.filepath,header = TRUE, stringsAsFactors = FALSE)
  baVermis.dat <- read.csv(baVermis_62.filepath,header = TRUE, stringsAsFactors = FALSE)
  sample.dat <- read.csv(samplefilePath,header = TRUE, stringsAsFactors = FALSE)
  
  ## read chromosome length file
  chrInfo = read.csv(file=chrFile, sep="\t", stringsAsFactors = FALSE)
  chrInd <- 1  
  while(chrInd <=25 ){
    #filter the datasets based on the chromosome ID
    chr.ba9_81.dat <- ba9_81.dat[which(ba9_81.dat[,1]==paste0("chr",chrInd)),]
    chr.ba41_66.dat <- ba41_66.dat[which(ba41_66.dat[,1]==paste0("chr",chrInd)),]
    chr.baVermis.dat <- baVermis.dat[which(baVermis.dat[,1]==paste0("chr",chrInd)),]
    nIter <- max(c(nrow(chr.ba9_81.dat),nrow(chr.ba41_66.dat),nrow(chr.baVermis.dat)))
  
    out <- NULL
    chrID <- toString(chrInfo$chrID[chrInd])
    chrLength <- chrInfo$chr_length[chrInd]
    
    ## Step-1: create the output dataframe and initialize with all 0's
    s <- seq(1, chrLength, binSize) 
    e <- s + 199
    df <- data.frame(matrix(0L, ncol = ncol(ba9_81.dat)+ncol(ba41_66.dat)+ncol(baVermis.dat)-6, nrow = length(s)))
    
    ## construct column header
    header <- c("chr","start","end", paste0(names(ba9_81.dat)[!names(ba9_81.dat) %in% c("chr","start","end")],".ba9"))
    header <- c(header, paste0(names(ba41_66.dat)[!names(ba41_66.dat) %in% c("chr","start","end")],".ba41"))
    header <- c(header, paste0(names(baVermis.dat)[!names(baVermis.dat) %in% c("chr","start","end")],".baVermis"))
    names(df) <- header ## getting error
    
    ## fill the windows
    df$chr = paste0("chr",chrInd)
    df$start = s
    df$end = e
    ## iteration
    for(i in 1:nIter){
      ######################        Process ba9
      # get chr region info
      if(!is.na(chr.ba9_81.dat[i,]$start)){
        rStart = chr.ba9_81.dat[i,]$start
        rEnd = chr.ba9_81.dat[i,]$end
        
        q1 <- rStart%/%binSize
        r1 <- rStart%%binSize
        q2 <- rEnd%/%binSize
        r2 <- rEnd%%binSize
        
        rowStartPos <- if(r1 <= binSize*overlapCutoff) (q1 + 1) else (q1 + 2)
        rowEndPos <- if((binSize-r2) <= binSize*overlapCutoff) (q2 + 1) else q2
        
        colStartPos <- 4
        colEndPos <- ncol(ba9_81.dat)
        df[rowStartPos:rowEndPos,colStartPos:colEndPos] <- 1  ## 1 indicates that there is a peak, i.e. a non-zero value (could be filled with actual value)
      }
      ######################        Process ba41
      # get chr region info
      if(!is.na(chr.ba41_66.dat[i,]$start)){
        rStart = chr.ba41_66.dat[i,]$start
        rEnd = chr.ba41_66.dat[i,]$end
        
        q1 <- rStart%/%binSize
        r1 <- rStart%%binSize
        q2 <- rEnd%/%binSize
        r2 <- rEnd%%binSize
        
        rowStartPos <- if(r1 <= binSize*overlapCutoff) (r1 + 1) else (r1 + 2)
        rowEndPos <- if((binSize-r2) <= binSize*overlapCutoff) (r2 + 1) else r2
        
        colStartPos <- ncol(ba9_81.dat) + 1
        colEndPos <- ncol(ba9_81.dat) + ncol(ba41_66.dat) - 3
        df[rowStartPos:rowEndPos,colStartPos:colEndPos] <- 1  ## 1 indicates that there is a peak, i.e. a non-zero value (could be filled with actual value)
      }
      
      ######################        Process baVermis
       # get chr region info
      if(!is.na(chr.baVermis.dat[i,]$start)){
        rStart = chr.baVermis.dat[i,]$start
        rEnd = chr.baVermis.dat[i,]$end
        
        q1 <- rStart%/%binSize
        r1 <- rStart%%binSize
        q2 <- rEnd%/%binSize
        r2 <- rEnd%%binSize
        
        rowStartPos <- if(r1 <= binSize*overlapCutoff) (r1 + 1) else (r1 + 2)
        rowEndPos <- if((binSize-r2) <= binSize*overlapCutoff) (r2 + 1) else r2
        
        colStartPos <- ncol(ba9_81.dat) + ncol(ba41_66.dat)-3+1
        colEndPos <- ncol(ba9_81.dat)+ncol(ba41_66.dat)+ncol(baVermis.dat)-6
        df[rowStartPos:rowEndPos,colStartPos:colEndPos] <- 1  ## 1 indicates that there is a peak, i.e. a non-zero value (could be filled with actual value)
      }
      
      if(i%%100 == 0) message(paste0("Iteration: ",i, " done"))
    }
    ## output the dataframe
    write.csv(df,file=paste0("df",chrInd,".csv"),row.names=F)
    message(paste0("chromosome ", chrInd, " has been completed!!"))
    chrInd <- chrInd + 1
  }
}
accetylationDat(ba9_81.filepath, ba41_66.filepath, baVermis_62.filepath, samplefilePath, chrFile, 200, overlapCutoff = 0.05, outputPath)
