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

### 1.1.6 Exporting transcript data for all the samples
- Use Rstudio and install 'ballgown' library.
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ballgown")
```

- Load all the abundance data (for more additional handy syntaxes see the 'ballgown' vignette)
```
## ----makebgobj, message=FALSE--------------------------------------------
library(methods)
library(ballgown)
library("data.table")
data_directory = "/Volumes/Seagate/STAR_Output/StringTieAbundance/"
# make the ballgown object:
bg = ballgown(dataDir=data_directory, meas='all', samplePattern="")
## ----get transcript spike-in (FPKM is the value we are interested in) ---
transcript_fpkm = texpr(bg, 'FPKM')
whole_tx_table = texpr(bg, 'all')
## get the transcript info, and output it
pref = whole_tx_table[,c(2,4,5)]
pref[,1]=paste0("chr",pref[,1])
newPref = cbind(pref,paste(pref[,1],pref[,2],pref[,3], sep="_"),".")
colnames(newPref) = c("chr","start","end","feature.id","strand")
fwrite(newPref,paste0(data_directory,"stringTie.Transcript.SpikeIns.bed"),col.names=F,quote=F,row.names=F)
## print the new StringTie SpikeIn (feature) values
newMat = cbind(pref,transcript_fpkm)
fwrite(newMat,paste0(data_directory,"stringTie.Transcript.SpikeIns.csv"),sep="\t",quote=F,row.names=F)
```
### 1.1.7 Bin the transcript abundance data

- "transcript_fpkm" or "transcript_cov" matrices contains transcript-level fpkm or coverage values in all the samples in all the chromosomes. Hence, we can apply "AccetylationDat3.r" (may need slight modification) code for binning these datasets into chromosome-wise files.  
- definition:
~cov: The average per-base coverage for the transcript or exon.
~FPKM: Fragments per kilobase of transcript per million read pairs. This is the number of pairs of reads aligning to this feature, normalized by the total number of fragments sequenced (in millions) and the length of the transcript (in kilobases).
~TPM: Transcripts per million. This is the number of transcripts from this particular gene normalized first by gene length, and then by sequencing depth (in millions) in the sample. A detailed explanation and a comparison of TPM and FPKM can be found here, and TPM was defined by B. Li and C. Dewey here.

- [update: 14/01/2019]: Use "transcript_fpkm" metrix, and convert it into a BED file for processing with the code "Accetylation_Bedtools.r"
