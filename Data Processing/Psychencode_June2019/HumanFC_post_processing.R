library(dplyr)
library(data.table)
# for HumanFC
# read the header (i.e. sample names)
setwd("/srv/scratch/z3526914/DeepBrain/Data/")
con <- file("mergedPeakHeightMatrix_HumanFC_filtered.bed","r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
dat <- fread("mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.bed", sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V293)) %>% select(-c(V293))
colnames(dat) <- header
fwrite(dat, file="mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.filtered.dat", sep="\t")
