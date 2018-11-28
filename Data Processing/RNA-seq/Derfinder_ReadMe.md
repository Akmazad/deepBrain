## 2 Derfinder: a R package for finding expressed region in RNA-seq data
### 2.1 Installation:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("derfinder", version = "3.8")
```
### 2.2 Usage: 
1) "windowMat.r" contains all the functions required to RNA-seq preparation (note: this code should be scaled up for running on raijin, including reading files from remote host()
2) "sample_code1.r" contains very basic codes demonstrating "derfinder" mainly
