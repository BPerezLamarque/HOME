example_great_apes_microbiota <-
function(name,path){ 
  data(great_apes_microbiota)
  write.tree(host_tree,file = paste(path,"/host_tree_",name,".tre",sep=""))
  write.dna(alignment_OTU892624276,paste(path,"/alignment_",name,"_OTU892624276.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(alignment_OTU892624276))
  write.dna(alignment_OTU47610657,paste(path,"/alignment_",name,"_OTU47610657.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(alignment_OTU47610657))
  write.dna(alignment_OTU733943228,paste(path,"/alignment_",name,"_OTU733943228.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=ncol(alignment_OTU733943228))
}
