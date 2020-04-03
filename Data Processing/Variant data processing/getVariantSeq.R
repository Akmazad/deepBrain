
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
# file-based processing: deepsea supple table 5 (noncoding GRASP eQTLs and negative variants sets)
filename="deepsea_supple_tabale5"
dat = as.data.frame(fread(paste0(filename, ".csv"),header = T))
dat[,8] = ifelse(dat[,8] == "eQTL",1,0)
dat = dat[,c(2,3,4,5,8)]
colnames(dat) <- c("Chr","Pos","Ref","Alt","Label")
dat$Strand = "+"

# invoke function for sequence retreival and save
fwrite(getVariabntSeq(dat, flankingLength), file = paste0(filename, "_fastaseq_refOnly.csv"), sep = ",")

