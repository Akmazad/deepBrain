# Input:  binsize, flankingLength, howManyBins (default: 1M)
#         mergedPeaks (mergedPeakHeightMatrix_HumanFC_filtered.bed, 
#         mergedPeakHeightMatrix_EpiMap_filtered.bed, 
#         Brain_CagePeaks_filtered.bed,
#         final.tf.bed)
# Host: RNA machine/KATANA 
#       RNA:    Data dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#               Script dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#       KATANA: Data dir: /srv/scratch/z3526914/DeepBrain/Data
#               Script dir: /srv/scratch/z3526914/DeepBrain/Scripts
options(scipen=999) # prevent auto-scientific notations of numeric values

library(dplyr)
library(data.table)
source = "HumanFC"
# ------------- Set Bin-flanking configuration
binSize = 400
flanking = 2*binSize     # can be arbitrarily given
baseDir = "/srv/scratch/z3526914/DeepBrain/Data/"
dataDir = paste0(baseDir, source, "/",binSize,"_",flanking,"/")
dir.create(dataDir, recursive=T) # create the directory if doesn't exists

# ------------- Read random Bins of given size
bin_outputFile = paste0("hg19_bins_", binSize,"bp","_rand")

# ---------------------------------------------- HumanFC
# --------------- IntersectBED with bins
humanFC_inputFile = "mergedPeakHeightMatrix_HumanFC_filtered"
# Note: "humanFC_inputFile" should exist within the "dataDir"
humanFC_outputFile = paste0(humanFC_inputFile,"_randBins.overlaps")
message("Intersect BED HumanFC:",appendLF=F)
system2('intersectBed', 
        paste('-wao -f 0.05 -a', paste0(baseDir,bin_outputFile,".bed"), '-b', paste0(baseDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# --------------- Drop unwanted fields
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".dropped")
message("Drop columns in HumanFC:",appendLF=F)
system2('cut', 
        paste('-f1-4,9-297',paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# --------------- Filter and aggregate (max) bins with same IDs
message("For multiple overlap with the same bin, pick the max-overlapped one: in HumanFC:",appendLF=F)
con <- file(paste0(baseDir,"mergedPeakHeightMatrix_HumanFC_filtered.bed"),"r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".filtered")
dat <- fread(paste0(dataDir,humanFC_inputFile,".bed"), sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V293)) %>% select(-c(V293))
colnames(dat) <- header
fwrite(dat, file=paste0(dataDir,humanFC_outputFile,".bed"), sep="\t")
message("Done",appendLF=T)

# --------------- Replace dots with 0
message("Replace dots with 0 (comes from intersectBed -wao): in HumanFC:",appendLF=F)
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".fixed")
system2('sed', 
        paste('\'s/\\./0/g\'',paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# --------------- Sort bins
message("Sort bins: in HumanFC:",appendLF=F)
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".sorted")
system2('sort', 
        paste('-k 1,1 -k2,2n',paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# --------------- non-Zero bins extraction in Shell (using awk): three file outputs
# non-Zero bins extraction in Shell (using awk): three file outputs
message("Getting Non-zero genomic bins:",appendLF=F)
# HumanFC
inputFile = "mergedPeakHeightMatrix_HumanFC_filtered_randBins.overlaps.dropped.filtered.fixed.sorted"
outputFile = "HumanFC_randBins_nonZero.binInfo"
system2('awk', 
        paste('-F \'\t\' \' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }\'',paste0(dataDir,inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# --------------- Extract genomic bins and lables and combine in a single file: one File output
# Extract genomic bins and lables and combine in a single file: one File output
library("BSgenome.Hsapiens.UCSC.hg19")
hg <- BSgenome.Hsapiens.UCSC.hg19
flankingLength <- flanking
inputFile <- outputFile
DNAoutputFile <- "HumanFC_randBins_nonZero.bin.Seq"
nonZerobins <- fread(paste0(dataDir,inputFile, ".bed"), sep="\t", header=F)
colnames(nonZerobins) <- c("chr","start","end","id")
fwrite(nonZerobins, file=paste0(dataDir,inputFile, ".bed"), sep="\t", quote=F, row.names=F)
seq <- getSeq(hg,nonZerobins$chr,start=nonZerobins$start-flankingLength,end=nonZerobins$end+flankingLength)
# seq <- tryCatch({
#   getSeq(hg,nonZerobins$chr,start=nonZerobins$start-flankingLength,end=nonZerobins$end+flankingLength)
# },error=function(e){
#   DNAString(paste(rep('N',2*flankingLength+binSize),collapse = ''))
# }
# )
seq <- as.character(as.data.frame(seq)[[1]])
nonZerobins.seq <- cbind(nonZerobins,seq)
colnames(nonZerobins.seq) <- c(colnames(nonZerobins),"dna.seq")
fwrite(nonZerobins.seq, file=paste0(dataDir, DNAoutputFile,".bed"), sep="\t", row.names=F, quote=F)

# --------------- Extract labels for non-zero bins
message("Extract labels for non-zero bins from each data files (HumanFC, EpiMap and ENCODE_TFs):",appendLF=F)
# HUMANFC
inputFile1 = "HumanFC_randBins_nonZero.binInfo.bed"
inputFile2 = "mergedPeakHeightMatrix_HumanFC_filtered_randBins.overlaps.dropped.filtered.fixed.sorted.bed"
HumanFC_nonzero_outputFile = "HumanFC_randBins_nonzero_labels"
system2('awk', 
        paste0('-F \"\t\" \'FILENAME==\"',
               paste0(dataDir,inputFile1),
               '\"{A[$1$2$3]=$1$2$3} FILENAME==\"',
               paste0(dataDir,inputFile2),
               '\"{if(A[$1$2$3]==$1$2$3){print}}\' ',
               paste0(dataDir,inputFile1),' ',
               paste0(dataDir,inputFile2)), 
        stdout=paste0(dataDir,HumanFC_nonzero_outputFile,".bed"), 
        wait=T)

message("Done",appendLF=T)

# HumanFC_nonzero_outputFile = "HumanFC_randBins_nonzero_labels"

# --------------- Print all the nonzero bins with sequence and labels
message("Print all the nonzero bins with sequence and labels:",appendLF=F)
outputFile <- "HumanFC_randBins_nonZero.bin.Seq_Labels"
dna.dat <- fread(paste0(dataDir, DNAoutputFile, ".bed"), sep="\t", header=T)
human <- fread(paste0(dataDir, HumanFC_nonzero_outputFile, ".bed"), sep="\t", header=T)
human$id <- paste0(human$chr, "_", human$start, "_", human$end)

# output <- cbind(human, epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
output_full <- cbind(dna.dat[which(dna.dat$id %in% human$id), ], human[,-c(1:4)])
colnames(output_full) <- c(colnames(dna.dat), colnames(human)[-c(1:4)])
fwrite(output_full,file=paste0(dataDir, outputFile, ".bed"), sep="\t", row.names=F, quote=F)

message("Done",appendLF=T)

# --------------- split (train and test (default: 'chr'))
# split (train and test (default: 'chr')) and numpy in Shell (python code): four File outputs
# remember to explicitely add run permission to the python script (chmod +x split_and_Numpy_V2.py)
system2("python", 
        paste(paste0(dataDir,"split_and_Numpy_V2.py"), 
              "--datadir", dataDir,
              "--datafilename HumanFC_randBins_nonZero.bin.Seq_Labels",
              "--valid_chr_id 1", 
              sep = ' '), 
        wait=T)

# Train and Test with DL model (5 epochs) in Shell and output a single accuracy (median) of all labels from final epoch's vallidation
# Use Google CoLab
