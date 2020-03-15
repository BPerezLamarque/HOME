# Function to create adding new host tips with branch lengths close to 0 (if several OTU sequences per host species)


library(ape)
library(phytools)

# Load your host tree and one OTU alignment

host_tree <- read.tree("host_tree_name.tre") # host tree with tips corresponding to host species names
# (warning: a host species name should not be a subset of another host species name in the following code)

alignment <- read.dna("alignment_name_OTU1.fas", as.character = T, format="fasta") # fasta alignment for OTU1 with faster header corresponding to host species name, plus an sequence identifier (if several OTU sequences are present per host species)



# This function creates new tips corresponding to sequences in the OTU alignment: it creates as many tips per host species in the tree as there are OTU sequences corresponding to each species. The new tips are created with branch lengths close to 0. 
add_host_tips <- function(host_tree, alignment){
  host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length)
  tip_labels <- host_tree$tip.label[order(nchar(host_tree$tip.label), decreasing = TRUE)]
  list_reads <- c()
  for (i in 1:length(tip_labels)){
    reads <- grep(tip_labels[i], rownames(alignment))
    reads <- reads[!reads %in% list_reads]
    list_reads <- c(list_reads, reads)
    if (length(reads)>1){
      if (!tip_labels[i] %in% rownames(alignment)){ host_tree$tip.label[which(host_tree$tip.label==tip_labels[i])] <-  rownames(alignment)[reads[1]]
      reference_tip <- rownames(alignment)[reads[1]] }else{ reference_tip <- tip_labels[i] }
      reads <- reads[which(!rownames(alignment)[reads] %in% host_tree$tip.label)]
      for (j in 1:length(reads)){
        host_tree <- bind.tip(host_tree, tip.label=rownames(alignment)[reads[j]], edge.length=NULL, where=which(host_tree$tip.label==reference_tip), position=0.001)
      }
    }
  }
  host_tree <- drop.tip(host_tree, tip = host_tree$tip.label[!host_tree$tip.label %in% rownames(alignment)])
  host_tree$edge.length[host_tree$edge.length==0] <- 0.001
  return(force.ultrametric(host_tree))
}

provided_tree <- add_host_tips(host_tree, alignment)

# the provided_tree can directly be used to run HOME. 

