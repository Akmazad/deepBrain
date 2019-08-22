# Data preprocessing for Deep learning training by analysing RNA-seq, CHIP-seq, ATAC-seq, and Transcription factor peaks (also CHIP-seq from UCSC-DCC)
Our pipeline considers only those chromosomal bins for DL training that has at least one TF features activated. Hence, we first constructed bins with at least one TF signal found, and after that augmented other features from e.g. EpiMap, HumanFC, or other sources.
## 1.1 StringTie 
To get TFs that are expressed in brain tissue, we've analysed RNA-seq data ([```Sun et al.```](https://www.sciencedirect.com/science/article/pii/S0092867416314519)) using [```Stringtie```](https://ccb.jhu.edu/software/stringtie/). Its a unix-based tool for finding expressed transcriptomic, exonic, and intronic chromosomal regions in the RNA-seq data.

### 1.1.1 Installation:
```
git clone https://github.com/gpertea/stringtie
cd stringtie
make release 
```

### 1.1.2 Assemble a batch of sample (BAM) files 
- To assemble the mapped reads into transcripts from all the sample (BAM) files, create a bash script (see "StringTieAssembly.sh") listing commands for all the samples and the output file paths. An example command for a single BAM file looks like following:
```
stringtie /Volumes/Seagate/STAR_Output/AN00493_ba41_42_22/AN00493_ba41_42_22Aligned.out.sorted.bam -l AN00493_ba41_42_22 -p 8 -G /Users/rna/Homo_sapiens.GRCh38.94.chr.gtf -o /Volumes/Seagate/STAR_Output/stringtieAssembly/AN00493_ba41_42_22.gtf
```

- To help running this scripts, paste this into home directory, for example, and add that location to the PATH variable using this command: [```$ export PATH=$PATH":$YourHomeDirectory"```]. 
- Then give execute permission to that script by [```$ chmod +x StringTieAssembly.sh```], which will enable running StringTie as a command.
- Reload the bash profile [```$ ~/.bash_profile```]
- Then run that script with [```./StringTieAssembly.sh```] to get the assembled files in .gtf format
- This takes help from the annotation file downloaded from here (```ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz```)

### 1.1.3 Merge all transcripts from the different samples
"StringTieMergeList.txt" contains all the output file names (.gtf) files
```
stringtie --merge -p 8 -G Homo_sapiens.GRCh38.94.chr.gtf -o /Volumes/Seagate/STAR_Output/stringtie_merged.gtf /Volumes/Seagate/STAR_Output/StringTieMergeList.txt
```

### 1.1.4 Count how many transcripts? [off-the-pipeline]
```
cat /Volumes/Seagate/STAR_Output/stringtie_merged.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l
```

### 1.1.5 Estimate transcript abundance [off-the-pipeline]
- Run "StringTieAbundance.sh". An examle line is as follows:
```
stringtie -e -B -p 8 -G /Volumes/Seagate/STAR_Output/stringtie_merged.gtf -o /Volumes/Seagate/STAR_Output/StringTieAbundance/AN00493_ba41_42_22/AN00493_ba41_42_22.gtf /Volumes/Seagate/STAR_Output/AN00493_ba41_42_22/AN00493_ba41_42_22Aligned.out.sorted.bam
```
- This will create a separate folder for each samples, where each folder will contain: 'i_data.ctab' (intron data), 'e_data.ctab' (exon data), 't_data.ctab' (transcript data), and their corresponding indices.

### 1.1.6.1 Exporting transcript data for all the samples
- Use Rstudio and install 'ballgown' library.
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ballgown")
```

- Load all the abundance data (for more additional handy syntaxes see the 'ballgown' vignette)
```r
## ----makebgobj, message=FALSE--------------------------------------------
library(methods)
library(ballgown)
library("data.table")
data_directory = "/Volumes/Seagate/STAR_Output/StringTieAbundance/"
# make the ballgown object:
bg = ballgown(dataDir=data_directory, meas='all', samplePattern="")
## ----get transcript spike-in (FPKM is the value we are interested in) ---
transcript_fpkm = texpr(bg, 'FPKM')
```
# 2: Select TFs using transcripts from Stringtie RNA-seq analysis
### 2.1 Determine a threshold for the percentage of samples with non-zero RNA-seq samples [Side analysis]
```r
fpkm_val_th=0.0 
transcript_fpkm = texpr(bg, 'FPKM')
whole_tx_table = texpr(bg, 'all')
plot_dat<-NULL
for(fpkm_perc_th in seq(0.1,1,0.05)){
    filtered.row=which(rowMeans(transcript_fpkm > fpkm_val_th) >= fpkm_perc_th)
    print(length(filtered.row))
    pref = whole_tx_table[filtered.row,c(2,4,5)]
    pref[,1]=paste0("chr",pref[,1])
    pref = cbind(pref, whole_tx_table[filtered.row,c(3,10)])
    gene_names=pref$gene_name
    nCommonGenes = length(intersect(tf_genes,gene_names))
    plot_dat = rbind(plot_dat,cbind(fpkm_perc_th,nCommonGenes))
}
colnames(plot_dat)=c("fpkm_perc_th","nTFs")
plot_dat=as.data.frame(plot_dat)
pdf("rplot.pdf") 
plot(plot_dat$fpkm_perc_th,plot_dat$nTFs))
dev.off() 
```
### 2.2 Filter the transcription matrix (with val_th & perc_th)
```r
fpkm_val_th=0.0 # decided this threshold since TFs expression are very low in general
fpkm_perc_th=0.5    # decided after "rplot.pdf" still picks as many TFs (150/161) as found in UCSC

filtered.row=which(rowMeans(transcript_fpkm > fpkm_val_th) >= fpkm_perc_th)
transcript_fpkm=transcript_fpkm[filtered.row,]
whole_tx_table = texpr(bg, 'all')
## get the transcript info, and output it
pref = whole_tx_table[filtered.row,c(2,4,5)]
pref[,1]=paste0("chr",pref[,1])
pref = cbind(pref, whole_tx_table[filtered.row,c(3,10)])
colnames(pref) = c("chr","start","end","strand","gene_name")
newMat = cbind(pref,transcript_fpkm)
output_directory="/Volumes/Data1/PROJECTS/DeepLearning/Test/"
fwrite(pref,paste0(output_directory,"stringTie.Transcript.SpikeIns_filtered.bed"),sep="\t",quote=F,row.names=F)
# Draw the frequency distribution of the filtered expression matrix
pdf("FD_filtered_transcript_matrix.pdf") 
hist(transcript_fpkm)
dev.off() 
# save the expressed TF_genes : should be 128 TF genes
tf.genes.expr = intersect(tf_genes,pref$gene_name)
```
### 2.3 Download and extract all the TF peak profiles (BED format; total 691) from the UCSC DCC portal [one-off]
```r
library(RCurl)
library(XML)
url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/"
doc <- htmlParse(url)
ucsc.tf.profileName <- as.data.frame(xpathSApply(doc, "//a/@href"))
ucsc.tf.profileName <- as.data.frame(ucsc.tf.profileName[-c(1:7),]) ## removing initial junk links
colnames(ucsc.tf.profileName)="filenames"
allFilePaths=paste0(url,ucsc.tf.profileName$filenames)
encodeFileOutDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCFiles/"

## download files : iterate for all, using 'apply' function
download.Save.file <- function(url,f) {
    tmp <- tempfile()
    download.file(f,tmp)
    data <- read.table(gzfile(tmp), sep="\t", header=F, stringsAsFactors=FALSE)
    #fName <- paste0(substr(f, nchar(url)+1, nchar(f)-2),"bed")
    fName <- substr(f, nchar(url)+1, nchar(f)-3)
    write.table(data,paste0(encodeFileOutDir,fName), quote = F, row.names = F, col.names = F, sep = "\t")
    unlink(tmp)
}
apply(as.array(allFilePaths), MARGIN = 1, FUN = function(x) download.Save.file(url,x))
```

### 2.4 Select peak profiles for expressed TFs (should be 595/691) 
```r
lookupDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/"
ucsc.tf.profileName=read.csv(paste0(lookupDir,"UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv"),header=T)
expr.ucsc.tf.profileName=ucsc.tf.profileName[ucsc.tf.profileName$Factor %in% tf.genes.expr,]$fileTitle
expr.ucsc.tf.profileName=paste0(expr.ucsc.tf.profileName,".narrowPeak")
oldDir=paste0(lookupDir,"EncodeDCCFiles/")
newDir=paste0(lookupDir,"EncodeDCCExprMatchFiles/")
# copy the files to the new folder
file.copy(paste0(oldDir,expr.ucsc.tf.profileName), newDir)

##write.csv(expr.ucsc.tf.profileName,paste0(lookupDir,"expr.ucsc.tf.profileName.csv"),row.names=F,quote=F)
#processAfiles <- function(aFile,lookupDir){
#    tfProfile=read.table(paste0(lookupDir,"EncodeDCCFiles/",aFile),header=T,sep="\t")
#    write.table(tfProfile,paste0(lookupDir,"EncodeDCCExprMatchFiles/",aFile),quote=F,row.names=F,sep="\t")
#}
#apply(as.array(expr.ucsc.tf.profileName), MARGIN = 1, FUN = function(x) processAfiles(x,lookupDir))
##library(stringr)
##expr.ucsc.tf.profileName=apply(as.array(tf.genes.expr),MARGIN = 1, FUN = function(x) str_detect(ucsc.tf.profileName,str_to_title(x)))
```

### 2.5 Merge all TF peak profiles into a single matrix
This section is influenced by [```peak_processing_new.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_processing_new.R) code with necessary modifications.
#### 2.5.1 For each sample, add a column with the sample name and one with peakCoordinate_sample
```r
setwd('/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCExprMatchFiles/')
files <- list.files(getwd(),pattern="*.narrowPeak$")
for (j in c(1:length(files)))
{
  data=read.table(as.character(files[j]), sep="\t", header=FALSE)
  sample=gsub(".narrowPeak", "", files[j])
  data=cbind(data,sample,paste(data[,1],data[,2],data[,3], sample, sep="_"))
  write.table(data, paste0(files[j],".sampleID.bed"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  print(j)
}
```
#### 2.5.2 Merge peaks using BEDTOOLS
```sh
cd /Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCExprMatchFiles/
outdir=/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCC_PeakProccessedData/
# concatenate all peak files in one; runtime 3min
cat *.sampleID.bed >> $outdir/allPeaks.sampleID.bed
# sort and merge peaks (mergebed requires sorted input)
sortbed -i $outdir/allPeaks.sampleID.bed > $outdir/allPeaks.sampleID.sorted.bed
mergebed -c 12 -o collapse -i $outdir/allPeaks.sampleID.sorted.bed > $outdir/allPeaks.sampleID.merged.bed
# wc -l $outdir/allPeaks.sampleID.merged.bed
#725276
```
#### 2.5.3 Generate a data matrix with peak height info for merged peaks
```r
library(tidyr)
library(dplyr)
library(data.table)
rm(list=ls())
setwd('/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCExprMatchFiles/')
merged=read.table("/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCC_PeakProccessedData/allPeaks.sampleID.merged.bed", sep="\t")
files <- list.files(getwd(),pattern="*.narrowPeak$")
# uncollapse merged peak table to one entry per row
merged$V4=as.character(merged$V4)
merged_uncollapsed= merged %>% unnest(V4 = strsplit(V4, ","))
dim(merged_uncollapsed)
#[1] 11881686        4
# remove duplicate rows. not sure why they exist
merged_uncollapsed=distinct(merged_uncollapsed)
dim(merged_uncollapsed)
#[1] 11795482        4

# split column 4 into initialPeakID and sample name and add a mergedPeakID
initialPeakInfo=transpose(strsplit(merged_uncollapsed$V4, "_"))
merged_uncollapsed=cbind(merged_uncollapsed, 
                         paste(merged_uncollapsed$V1, merged_uncollapsed$V2, merged_uncollapsed$V3, sep="_"),
                         paste(initialPeakInfo[[1]], initialPeakInfo[[2]], initialPeakInfo[[3]], sep="_"),
                         initialPeakInfo[[4]])
colnames(merged_uncollapsed)=c("chr", "start", "end", "initialPeakID_Sample", "mergedPeakID", "initialPeakID", "Sample")
peaks=matrix(NA, nrow=length(unique(merged_uncollapsed$mergedPeakID)), ncol=length(files))
rownames(peaks)=unique(merged_uncollapsed$mergedPeakID)
colnames(peaks)=as.character(gsub(".narrowPeak", "", files))
# specify peakHeight column from the initial peak files (can differ between datasets), as well as initialPeak_SampleCol and sampleCol  (which were created in step 1).
peakHeightCol=7
sampleCol=11  
initialPeak_SampleCol=12  
for (j in c(1:length(files)))
{ # read in peak height data
  data=read.table(paste0(as.character(files[j]), ".sampleID.bed"), sep="\t", header=FALSE)
  # subset merged peaks by those with contributing peaks from sample j
  use_peaks=merged_uncollapsed[which(merged_uncollapsed$Sample%in%data[,sampleCol]),]
  use_peaks$peakHeight=data[match(use_peaks$initialPeakID_Sample, data[,initialPeak_SampleCol]), peakHeightCol]
  # sum all peaks from sampleJ by merged peak coordinates
  peak_sum=aggregate(use_peaks$peakHeight, by=as.data.frame(use_peaks$mergedPeakID), FUN="sum")
  # add data to the matrix            
  s=which(colnames(peaks)%in% as.character(data[, sampleCol]))
  p=match(rownames(peaks), peak_sum[,1])
  peaks[,s]=peak_sum[p,2]
  print(j)
}
save(peaks,file= "/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCC_PeakProccessedData/mergedPeakHeightMatrix_ENCODE_DCC_TFs.rda")

# Peak filtering (from peak_filtering.R)
outdir <- "/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCC_PeakProccessedData/"
filename <- "mergedPeakHeightMatrix_EncodeDCC_TFs"
val_th <- 0 # can be decided later
new_peaks <- ifelse(peaks>val_th, 1, 0)
rm(peaks)
new_peaks <- ifelse(!is.na(new_peaks), 1, 0) # replacing NAs with 0
binInfo <- as.data.frame(do.call(rbind,strsplit(rownames(new_peaks), "_")))
final.dat <- cbind(binInfo, rownames(new_peaks), new_peaks)
colnames(final.dat) <- c("chr","start","end", "id", colnames(new_peaks))
fwrite(final.dat, paste0(outdir, filename, "_binarized.BED"), row.names = F, quote = F, sep = "\t")
```
#### 2.5.4 convert it TF gene-symbol based matrix
```r
setwd('/Volumes/Data1/PROJECTS/DeepLearning/Test/')
tfProfFile = "UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv"
tfProf <- read.csv(tfProfFile, header = T)[,c(2,6)]
tfs = read.csv("tf.genes.expr_0.5_perc.csv", header = T) 
processAtf <- function(aTF,dat,tfProf){
  print(aTF)
  aTF.Prof = as.vector(tfProf[as.character(tfProf$Factor)==aTF,2])
  aTF.Prof.dat = as.data.frame(dat[,aTF.Prof])
  agg.val <- data.frame(x1 = rowSums2(as.matrix(aTF.Prof.dat)))
  colnames(agg.val)=aTF
  agg.val= as.data.frame(ifelse(agg.val>0,1,0))
  return(agg.val)
}
temp <- apply(as.matrix(tfs),1,FUN=function(x) processAtf(x,final.dat,tfProf)) %>% as.data.frame()
final.dat.tf <- cbind(final.dat[,1:4], temp) 
colnames(final.dat.tf) <- c("chr", "start", "end", "id", colnames(temp))
fwrite(final.dat.tf,file = "final.dat.tf", sep="\t", row.names=F, quote=F)

```
