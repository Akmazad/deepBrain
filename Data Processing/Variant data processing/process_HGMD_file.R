library(data.table)
library(dplyr)
library(tidyr)

complement <- function(X){
  if(X == 'A')
    return('T')
  else if(X == 'T')
    return('A')
  else if(X == 'C')
    return('G')
  else if(X == 'G')
    return('C')
}

setwd("C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\")
dat <- fread("HGMD_Search_Results_PromoterActivity.txt") %>% as.data.frame() %>% dplyr::select(6,11)
dat1 <- dat %>% tidyr::separate(1, c("SNP_chr","SNP_start","strand"), ":") %>% na.omit()
dat1$SNP_end <- dat1$SNP_start
dat2 <- (dat1 %>% tidyr::separate(4, c("X1","X2"), "\\["))[,-4] %>%
  tidyr::separate(4, c("X3","X4"), "\\]") %>% dplyr::select(-5) %>% 
  tidyr::separate(4, c("Ref","Alt"), "\\/")
change <- which(dat2$strand == '-')
dat3 = dat2
dat3[change,4] <- sapply(dat3[change,4], function(x){
  return(complement(x))
})
dat3[change,5] <- sapply(dat3[change,5], function(x){
  return(complement(x))
})
dat3$SNP_id = paste0(dat3$SNP_chr, "_", dat3$SNP_start,  "_", dat3$SNP_end)
dat3 = dat3[,c("SNP_chr","SNP_start","SNP_end","SNP_id","Ref","Alt", "strand")]
fwrite(dat3, file = "HGMD_Search_Results_PromoterActivity_processed.csv", sep = "\t", col.names = F)
