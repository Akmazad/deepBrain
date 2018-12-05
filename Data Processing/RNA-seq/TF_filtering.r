dat <- read.csv("/Volumes/Seagate/STAR_Output/StringTie_all_Transcript_info.csv")
## download TF gene symbols for human from "http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.xlsx"
tf.dat <- read.csv("/Volumes/Seagate/STAR_Output/TFCheckpoint_download_180515.csv")
tf.human.dat <- tf.dat[which(tf.dat$entrez_human != 0),]
dbtf.human.dat <- tf.human.dat[which(tf.human.dat$DbTF == 'yes'),]

## filter human tfs only
dat <- dat[which(dat$gene_name %in% dbtf.human.dat$gene_symbol),]
write.csv(dat, file="/Volumes/Seagate/STAR_Output/StringTie_humanTFsOnly.csv",row.names=F)
