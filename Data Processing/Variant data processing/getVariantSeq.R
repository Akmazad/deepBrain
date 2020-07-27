getVariabntSeq <- function(dat, flankingLength){
  require(data.table)
  require(BSgenome.Hsapiens.UCSC.hg19)
  require(dplyr)
  require(stringr)
  
  hg <- BSgenome.Hsapiens.UCSC.hg19
  start = dat$Pos - flankingLength
  end = dat$Pos + (flankingLength - 1)
  # chrName= paste0("chr",as.character(dat$Chr))
  
  # ref seq
  fasta.seq.ref = BSgenome::getSeq(hg, dat$Chr, start=start, end=end, strand = as.character(dat$Strand))
  fasta.seq.ref = as.data.frame(fasta.seq.ref)[,1]
  # print(substring(fasta.seq.ref, (flankingLength + 1), (flankingLength + 1)))\
  message("ref Fasta seq retrieval: [DONE]")

  # variant seq
  fasta.seq.var = fasta.seq.ref
  substr(fasta.seq.var, (flankingLength + 1), (flankingLength + 1)) <- as.character(dat$Alt)
  message("variant Fasta seq retrieval: [DONE]")
  
  # print(substring(fasta.seq.var, (flankingLength + 1),(flankingLength + 1)))
  # substring(fasta.seq,(flankingLength + 1),(flankingLength + 1))

  # dat = cbind(fasta.seq.ref, fasta.seq.var, dat$Label)
  # colnames(dat) <- c("refDNAseq", "varDNAseq", "Label")
  dat = cbind(dat$Ref, dat$Alt, flankingLength, fasta.seq.ref, dat$Label)
  colnames(dat) <- c("Ref","Alt", "Pos", "refDNAseq", "Label")
  return(dat)
}

# start
flankingLength = 500
library(data.table)
library(dplyr)

dir = "/srv/scratch/z3526914/DeepBrain/Data/"
# dir = "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\"

# True Positives (i.e. Functional SNPs - Brain-based)
filename="brain_specific_eQTL"
dat = fread(paste0(dir, filename, ".txt")) %>% as.data.frame()
dat = dat[,c(19,11,20,21,6)]
colnames(dat) <- c("Chr","Pos","Ref","Alt","Strand")
# find unique snps
dat = unique(dat)
fwrite(dat, paste0(dir,"SNPS.BED"), header=F)   # this list has both coding and non-coding SNPs)
message("Intersect BED to get non-coding based SNPs:",appendLF=F)
system2('intersectBed', 
        paste('-V -A ', paste0(dir,"SNPS.BED"), '-B', paste0(dir,"CODING.EXONS.BED"), sep=' '), 
        stdout=paste0(dir,"NONCODING.SNPS.BED"), 
        wait=T)
message("Done",appendLF=T)

dat$Lable = "1" # all TP
fwrite(getVariabntSeq(dat, flankingLength), file = paste0(dir, filename, "_fastaseq.csv"), sep = ",")

# True Negative (i.e. non-Functional SNPs - from deepSea supple table 5)
filename="deepsea_supple_tabale5"
# anything other than "eQTL" is negative
dat = as.data.frame(fread(paste0(filename, ".csv"),header = T)) %>% dplyr::filter(label != "eQTL")
# dat[,8] = ifelse(dat[,8] == "eQTL",1,0)
dat = dat[,c(2,3,4,5,8)]
colnames(dat) <- c("Chr","Pos","Ref","Alt","Label")
dat$Strand = "+"

# invoke function for sequence retreival and save
fwrite(getVariabntSeq(dat, flankingLength), file = paste0(dir, filename, "_fastaseq_Negatives.csv"), sep = ",")
