#########   SELECT THE CORE OTU FOR RUNNING HOME #########   

rm(list=ls())

setwd("YOUR_WORKING_DIRECTORY")

library(ape)

table <- read.table("OTU_table_16S_1f_filtered.txt",header=T)

colnames(table)[1:11] <- c("abundance","cloud","OTU","length","abundance_2","chimera","spread","sequence","identity","reference","taxonomy")

sum(table[,12:ncol(table)]) 

### YOU SHOULD HERE REMOVED THE OTU CONTAMINANTS 

# If an OTU has less than 5 sequences in a sample, it is not considered (i.e., it's considered as mising in the sample)
for (i in 12:ncol(table)){table[which(table[,i]<5),i] <- 0}


threshold_core=10 ### Number of samples in which the OTU is found to be considered as a core OTU



length(which(table$spread>=threshold_core))
table$taxonomy[which(table$spread>=threshold_core)]

table_core <- table
table_core <- table_core[which(table_core$spread>=threshold_core),]

hist(table_core$spread)

write.table(table_core,file=paste0("OTU_table_core_OTUs.csv"),sep=";",quote=F,row.names = F)


for (i in 1:nrow(table_core)){
  list <- colnames(table_core)[11+which(table_core[i,12:ncol(table_clean)]>=5)] 
  write.table(list,file=paste0("core_OTU/list_sample_OTU_",table_core$OTU[i],".txt"),sep="\n",quote=F,row.names = F,col.names = F)
}

write.table(table_core$OTU,file=paste0("core_OTU/list_core_OTU.txt"),sep="\n",quote=F,row.names = F,col.names = F)



