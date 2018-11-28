## 1: Preparing RNA-seq data from Brain tissue for Deep learning applications
# 1.1 Derfinder
1) "windowMat.r" contains all the functions required to RNA-seq preparation (note: this code should be scaled up for running on raijin, including reading files from remote host()
2) "sample_code1.r" contains very basic codes demonstrating "derfinder" mainly

# 1.2 StringTie: a unix-based tool 

# Installation:
```
git clone https://github.com/gpertea/stringtie
cd stringtie
make release 
```

# Run a batch of sample (BAM) files 
- To assemble the mapped reads into transcripts from all the sample (BAM) files, create a bash script (see "StringTieAssembly.sh") listing commands for all the samples and the output file paths. 
- To help running this scripts, paste this into home directory, for example, and add that location to the PATH variable using this command: [$ export PATH=$PATH":$YourHomeDirectory"]. 
- Then give execute permission to that script by [$ chmod +x StringTieAssembly.sh], which will enable running StringTie as a command.
- Then run that script with [./StringTieAssembly.sh] to get the assembled files in .gtf format
- This takes help from the annotation file downloaded from here (ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz)

