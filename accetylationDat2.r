## preconditions: 
## 1. A window (i.e. [wStart:wEnd]) will not completely engulf the chromosome region (i.e. [rStart:rEnd])
## 2. Hence, there WILL NOT be more than 1 overlap for each windown
## 3. overlapCutoff value is between [0,1] and applied on the chromosomal region (i.e. [rStart:rEnd])

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

accetylationDat <- function(ba9_81.filepath, ba41_66.filepath, baVermis_62.filepath, samplefilePath, chrFile, binSize, overlapCutoff, outputPath){
  library("reproducible")
  ### load files
  ba9_81.dat <- read.csv(ba9_81.filepath,header = TRUE, stringsAsFactors = FALSE)
  ba41_66.dat <- read.csv(ba41_66.filepath,header = TRUE, stringsAsFactors = FALSE)
  baVermis.dat <- read.csv(baVermis_62.filepath,header = TRUE, stringsAsFactors = FALSE)
  sample.dat <- read.csv(samplefilePath,header = TRUE, stringsAsFactors = FALSE)
  
  ## read chromosome length file
  chrInfo = read.csv(file=chrFile, sep="\t", stringsAsFactors = FALSE)
  chrInd <- 1
  out <- NULL
  
  while(chrInd <=1 ){
    chrID <- toString(chrInfo$chrID[chrInd])
    chrLength <- chrInfo$chr_length[chrInd]
    wStart <- 18910
    binID <- 1
    #while (wStart < chrLength){  
    while (wStart < 19110){
      wEnd <- if(wStart + binSize - 1 <= chrLength) (wStart + binSize - 1) else chrLength
      chrName <- as.character(lapply(strsplit(dat[,1],"_"), `[[`, 1))
      
      #######################    FOR "ba9_81" brain region    #######################
      dat <- ba9_81.dat
      rStart <- dat[,2] - dat[,2]*overlapCutoff
      rEnd <- dat[,3] + dat[,3]*overlapCutoff
      cond1 <- "&"(chrName==paste0("chr",chrInd), "&"(dat[,2] <= wStart, dat[,3] >= wEnd)) ## whole window contained
      cond2 <- "&"(chrName==paste0("chr",chrInd), "&"(rStart <= wStart, dat[,3] >= wEnd))   ## window start is before the chr start
      cond3 <- "&"(chrName==paste0("chr",chrInd), "&"(dat[,2] <= wStart, rEnd >= wEnd))   ## window start is before the chr start   ## window end is after the chr end
      cond <- "|"(cond1,cond2)
      cond <- "|"(cond,cond3)
      temp.ba9_81 <- ba9_81.dat[which(cond),]
      
      # if nothing found, make an empty row with zeros
      if(nrow(temp.ba9_81) < 1){  
        xx1 <- rep(0,length(temp.ba9_81))
        nam <- names(temp.ba9_81)
        temp.ba9_81 <- rbind(temp.ba9_81,xx1)
        names(temp.ba9_81) <- nam
      }
      # trim off some columns
      exc <- names(temp.ba9_81) %in% c("chr","start","end")
      temp.ba9_81 <- temp.ba9_81[!exc]
      names(temp.ba9_81) <- paste0(names(temp.ba9_81),".","ba9")
      
      #######################    FOR "ba41_66" brain region    #######################
      dat <- ba41_66.dat            
      rStart <- dat[,2] - dat[,2]*overlapCutoff
      rEnd <- dat[,3] + dat[,3]*overlapCutoff
      cond1 <- "&"(chrName==paste0("chr",chrInd), "&"(dat[,2] <= wStart, dat[,3] >= wEnd)) ## whole window contained
      cond2 <- "&"(chrName==paste0("chr",chrInd), "&"(rStart <= wStart, dat[,3] >= wEnd))   ## window start is before the chr start
      cond3 <- "&"(chrName==paste0("chr",chrInd), "&"(dat[,2] <= wStart, rEnd >= wEnd))   ## window start is before the chr start   ## window end is after the chr end
      cond <- "|"(cond1,cond2)
      cond <- "|"(cond,cond3)      
      temp.ba41_66 <- ba41_66.dat[which(cond),]
      
      if(nrow(temp.ba41_66) < 1){  # if nothing found, make an empty row with zeros
        xx1 <- rep(0,length(temp.ba41_66))
        nam <- names(temp.ba41_66)
        temp.ba41_66 <- rbind(temp.ba41_66,xx1)
        names(temp.ba41_66) <- nam
      }
      # trim off some columns
      exc <- names(temp.ba41_66) %in% c("chr","start","end")
      temp.ba41_66 <- temp.ba41_66[!exc]
      names(temp.ba41_66) <- paste0(names(temp.ba41_66),".","ba41-42-22")
      
      #######################    FOR "baVermis" brain region    #######################
      dat = baVermis.dat            
      rStart <- dat[,2] - dat[,2]*overlapCutoff
      rEnd <- dat[,3] + dat[,3]*overlapCutoff
      cond1 <- "&"(startsWith(dat[,1],paste0("chr",chrInd)), "&"(dat[,2] <= wStart, dat[,3] >= wEnd)) ## whole window contained
      cond2 <- "&"(startsWith(dat[,1],paste0("chr",chrInd)), "&"(rStart <= wStart, dat[,3] >= wEnd))   ## window start is before the chr start
      cond3 <- "&"(startsWith(dat[,1],paste0("chr",chrInd)), "&"(dat[,2] <= wStart, rEnd >= wEnd))   ## window start is before the chr start   ## window end is after the chr end
      cond <- "|"(cond1,cond2)
      cond <- "|"(cond,cond3)      
      temp.baVermis <- baVermis.dat[which(cond),]
      
      if(nrow(temp.baVermis) < 1){  # if nothing found, make an empty row with zeros
        xx1 <- rep(0,length(temp.baVermis))
        nam <- names(temp.baVermis)
        temp.baVermis <- rbind(temp.baVermis,xx1)
        names(temp.baVermis) <- nam      
      }
      # trim off some columns
      exc <- names(temp.baVermis) %in% c("chr","start","end")
      temp.baVermis <- temp.baVermis[!exc]
      names(temp.baVermis) <- paste0(names(temp.baVermis),".","vermis")
      
      ## See Precondition 1 and 2
      if(nrow(temp.ba9_81) > 1 || nrow(temp.ba41_66) > 1 || nrow(temp.baVermis) > 1) message("Should not happen !!")      
      
      aRow <- cbind(chrInd,wStart,wEnd,temp.ba9_81,temp.ba41_66)
      aRow <- cbind(aRow,temp.baVermis)
      
      if(!all(aRow == 0))
        out <- rbind(out,aRow)
      wStart <- wEnd + 1
      binID <- binID + 1
    }
    chrInd <- chrInd + 1
  }
  
  #write.csv(out,file = paste(outputPath, "\\outputAccetylation.csv"))
  write.csv(out,file = "outputAccetylation.csv", row.names=FALSE)
}

accetylationDat(ba9_81.filepath, ba41_66.filepath, baVermis_62.filepath, samplefilePath, chrFile, 200, overlapCutoff = 0, outputPath)
