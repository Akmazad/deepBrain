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
### 2.5 Intersect Bins with selected TF peaks
```r
############## Overlap Bins with TFB locations, with a min of 5% overlap ; done in shell using bedTools (can be embeded in R)
  # Step-1: create a shell script namely "tfbs_bins_overlap.sh" (see attached) within the "workingDir"
  # Step-2: register that script for running with appropriate permission under UNIX using "chmod u+x tfbs_bins_overlap.sh"
  # Step-3: Put following commands for Bedtools in that shell script which assumes the arguments should be passed from R
  
  ##### "tfbs_bins_overlap.sh" ########
  ## #!/bin/bash
  ## bedDir=$1
  ## fileDir=$2
  ## bins=$3/$4
  ## overlapCutof=$(echo $5| bc)
  ## for tfFile in "$fileDir"/*.narrowPeak
  ## do
    ## $bedDir/intersectBed -u -f $overlapCutof -a $bins.bed -b $tfFile > $tfFile.overlaps.bed
    ## filename="${tfFile##*/}"
    ## echo "$filename: [DONE]"
  ## done

  
  # Step-4: use system2 R function to run this script with arguments passed for the shell script
  bedDir="/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
  fileDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCExprMatchFiles"
  binFileDir="/Volumes/Data1/PROJECTS/DeepLearning/Test"
  binFile="hg19_bins_200bp"
  overlapCutoff=0.05
  
  message(paste0("Overlapping tf locations with Bin locations, with a min of ",overlapCutoff*100, "% overlap: "),appendLF=F)
  system2("./tfbs_bins_overlap.sh",
            paste(bedDir, 
            fileDir, 
            binFileDir, 
            binFile, 
            overlapCutoff,
            sep=" "))
  # this will create Three overlap bed files
 
 message("Done",appendLF=T)

```
### 2.6 Merge binned TF profiles with other features
