## Example run 

library(ape)
library(phytools)
library(HOME)

name="TERMITE"
name_index=c("K02586")

alignment <- read.dna("alignment_TERMITE_K02586.fas", format = "fasta", as.character = T)
host_tree <- read.tree("host_tree_TERMITE.tre")

# Add tree tips with close to zero branch lengths to match number of microbial sequences
add_host_tips <- function(host_tree, alignment){
  host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length)
  tip_labels <- host_tree$tip.label[order(nchar(host_tree$tip.label), decreasing = TRUE)]
  list_reads <- c()
  for (i in 1:length(tip_labels)){
    reads <- grep(tip_labels[i], rownames(alignment))
    reads <- reads[!reads %in% list_reads]
    list_reads <- c(list_reads, reads)
    if (length(reads)>0){
      if (!tip_labels[i] %in% rownames(alignment)){ host_tree$tip.label[which(host_tree$tip.label==tip_labels[i])] <-  rownames(alignment)[reads[1]]
      reference_tip <- rownames(alignment)[reads[1]] }else{ reference_tip <- tip_labels[i] }
      reads <- reads[which(!rownames(alignment)[reads] %in% host_tree$tip.label)]
      if (length(reads)>0){
        for (j in 1:length(reads)){
          host_tree <- bind.tip(host_tree, tip.label=rownames(alignment)[reads[j]], edge.length=NULL, where=which(host_tree$tip.label==reference_tip), position=min(0.001,min(host_tree$edge.length)))
        }
      }
    }
  }
  host_tree <- drop.tip(host_tree, tip = host_tree$tip.label[!host_tree$tip.label %in% rownames(alignment)])
  host_tree$edge.length[host_tree$edge.length==0] <- 0.001
  return(force.ultrametric(host_tree,method = "extend"))
}

provided_tree <- add_host_tips(host_tree, alignment)


# Run HOME
HOME_model(name="TERMITE", name_index=c("K02586"), provided_tree = provided_tree, nb_tree=5000, lambda=seq(1,50), nb_random=10)


