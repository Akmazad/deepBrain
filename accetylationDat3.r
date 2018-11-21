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
  nR.ba9_81.dat <- nrow(ba9_81.dat)
  nR.ba41_66.dat <- nrow(ba41_66.dat)
  nR.baVermis.dat <- nrow(baVermis.dat)
  
  ## read chromosome length file
  chrInfo = read.csv(file=chrFile, sep="\t", stringsAsFactors = FALSE)
  chrInd <- 1  
  while(chrInd <=1 ){
    out <- NULL
    chrID <- toString(chrInfo$chrID[chrInd])
    chrLength <- chrInfo$chr_length[chrInd]
    
    ## Step-1: create the output dataframe and initialize with all 0's
    s <- seq(1, chrLength, binSize)
    df <- data.frame(matrix(0L, ncol = ncol(ba9_81.dat)+ncol(ba41_66.dat)+ncol(baVermis.dat)-6, nrow = length(s)))
    
  }
 
}
