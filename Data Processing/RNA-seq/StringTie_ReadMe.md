# 1: Preparing RNA-seq data from Brain tissue for Deep learning applications
## 1.1 StringTie 
A unix-based tool for finding expressed transcriptomic, exonic, and intronic chromosomal regions in the RNA-seq data

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

### 1.1.4 Count how many transcripts?
```
cat /Volumes/Seagate/STAR_Output/stringtie_merged.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l
```

### 1.1.5 Estimate transcript abundance
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
### 1.1.6.2 Load UCSC TF profile find TF genes that are expressed in the RNA-seq data
```r
## download Transcription Factor ChIP-seq Uniform Peaks from ENCODE/Analysis (source UCSC: http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeAwgTfbsUniform)
tf.dat <- read.csv("/Volumes/Data1/PROJECTS/DeepLearning/Test/UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv", header=F)
tf_genes <- unique(tf.dat[,2]) ## second column contains the gene-symbol
```
### 1.1.6.3 Determine a threshold for the percentage of samples with non-zero RNA-seq samples
```r
transcript_fpkm = texpr(bg, 'FPKM')
whole_tx_table = texpr(bg, 'all')
plot_dat<-NULL
fpkm_val_th=0.0 ## basicaly to ignore this parameter and reuse the existing filtering-code
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
### 1.1.6.4 Filter the transcription matrix (percentage threshold only)
```r
fpkm_val_th=0.0
fpkm_perc_th=0.5

filtered.row=which(rowMeans(transcript_fpkm > fpkm_val_th) >= fpkm_perc_th)
transcript_fpkm=transcript_fpkm[filtered.row,]
whole_tx_table = texpr(bg, 'all')
## get the transcript info, and output it
pref = whole_tx_table[filtered.row,c(2,4,5)]
pref[,1]=paste0("chr",pref[,1])
pref = cbind(pref, whole_tx_table[filtered.row,c(3,10)])
colnames(pref) = c("chr","start","end","strand","gene_name")
newMat = cbind(pref,transcript_fpkm)
output_directory="/Volumes/Data1/PROJECTS/DeepLearning/Test"
fwrite(pref,paste0(output_directory,"stringTie.Transcript.SpikeIns_filtered.bed"),sep="\t",quote=F,row.names=F)
# Draw the frequency distribution of the filtered expression matrix
pdf("FD_filtered_transcript_matrix.pdf") 
hist(transcript_fpkm)
dev.off() 
#save the TF_genes in the expression data
tf.genes.expr = intersect(tf_genes,pref$gene_name)
```
### 1.1.6.5 Overlap the filtered transctipts with the Bins
```r
bedDir="/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
fileDir="/Volumes/Data1/PROJECTS/DeepLearning/Test"
binFile="hg19_bins_200bp"
overlapCutoff=0.05
transctiptFile="stringTie.Transcript.SpikeIns_filtered"
message(paste0("Overlapping filtered transctipts with Bin locations, with a min of ",overlapCutoff*100, "% overlap: "),appendLF=F)
system2("./filteredTrascript_bins_overlap.sh",
            paste(bedDir, 
            fileDir, 
            binFile, 
            overlapCutoff,
            transctiptFile,
            sep=" "))
message("Done",appendLF=T)
## the output of this overlap is named as "stringTie.Transcript.SpikeIns_filtered.overlaps.bed", located in the "fileDir"
```
### 1.1.6.5.1 Further expression matrix (filtered) for only TFs that are expressed 
```r
tf_genes.rnaSeq =  newMat[which(newMat$gene_name %in% tf_genes),]
tf_genes.dat.ucscAcc = tf.dat[which(tf.dat$Factor %in% tf_genes.rnaSeq$gene_name),3]    # third column holds the UCSC accession number

```
### 1.1.6.5.2 Download and extract all the TF profiles from the UCSC DCC portal
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
    data <- read.table(gzfile(tmp), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    fName <- paste0(substr(f, nchar(url)+1, nchar(f)-2),"bed")
    write.table(data,paste0(encodeFileOutDir,fName), quote = F, row.names = F, sep = "\t")
    unlink(tmp)
}
apply(as.array(allFilePaths), MARGIN = 1, FUN = function(x) download.Save.file(url,x))
```
### 1.1.6.5.3 Pick the TF (expressed) profiles 
```r
lookupDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCFiles/"
ucsc.tf.profileName=list.files(path = lookupDir, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case =                                    FALSE, include.dirs = FALSE)
library(stringr)
expr.ucsc.tf.profileName=apply(as.array(tf.genes.expr),MARGIN = 1, FUN = function(x) str_detect(ucsc.tf.profileName,str_to_title(x)))
```
### 1.1.6b [test analysis] intersect with TFBS locations (Encode) 
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
  fileDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/UCSC_Encode_DCC"
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
### 1.1.6b [test analysis] intersect with Enhancer locations (PsychEncode) [IGNORED]
```r
############## Overlap StringTie features with Enhancer locations, with a min of 5% overlap ; done in shell using bedTools (can be embeded in R)
  # Step-1: create a shell script namely "rnaSeqEnhancer_Bed_ShellScript.sh" (see attached) within the "workingDir"
  # Step-2: register that script for running with appropriate permission under UNIX using "chmod u+x rnaSeqEnhancer_Bed_ShellScript.sh"
  # Step-3: Put following commands for Bedtools in that shell script which assumes the arguments should be passed from R
  
  ##### "rnaSeqEnhancer_Bed_ShellScript.sh" ########
  ## #!/bin/bash
  ## bedDir=$1	#first argument which is the Bedtools bin directory
  ## cd $2	#second argument which is the working directory
  ## rnaSeqBed=$3	#third argument which is the rnaSeq bed file name
  ## overlapCutof=$(echo $4| bc)	#fourth argument which is overlap cut-off
  ## enhancerBed=$5                 #fifth argument which is enhancerFileName
  ## $bedDir/intersectBed -u -f $overlapCutof -a $rnaSeqBed.bed -b $enhancerBed.bed > $rnaSeqEnhancer.overlaps.bed
  ## done
  
  # Step-4: use system2 R function to run this script with arguments passed for the shell script
  bedDir="/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
  workingDir="/Volumes/Data1/PROJECTS/DeepLearning/Test"
  strintTieBedFile="stringTie.Transcript.SpikeIns_full_binarized"
  overlapCutoff=0.05
  enhancerBedFile="DER-03_hg19_PEC_enhancers"
  message(paste0("Overlapping rnaSeq locations with enhancers locations, with a min of ",overlapCutoff*100, "% overlap: "),appendLF=F)
  system2("./rnaSeqEnhancer_Bed_ShellScript.sh",
            paste(bedDir, 
            workingDir, 
            strintTieBedFile, 
            overlapCutoff, 
            enhancerBedFile,
            sep=" "))
  # this will create Three overlap bed files
 
 message("Done",appendLF=T)
```

### 1.1.7 Bin the transcript abundance data

- "transcript_fpkm" or "transcript_cov" matrices contains transcript-level fpkm or coverage values in all the samples in all the chromosomes. Hence, we can apply "AccetylationDat3.r" (may need slight modification) code for binning these datasets into chromosome-wise files.  
- definition:
    - cov: The average per-base coverage for the transcript or exon.
    - FPKM: Fragments per kilobase of transcript per million read pairs. This is the number of pairs of reads aligning to this feature, normalized by the total number of fragments sequenced (in millions) and the length of the transcript (in kilobases).
    - TPM: Transcripts per million. This is the number of transcripts from this particular gene normalized first by gene length, and then by sequencing depth (in millions) in the sample. A detailed explanation and a comparison of TPM and FPKM can be found here, and TPM was defined by B. Li and C. Dewey here.

- [update: 14/01/2019]: Use "transcript_fpkm" metrix, and convert it into a BED file for processing with the code "Accetylation_Bedtools.r"
