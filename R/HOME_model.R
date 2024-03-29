HOME_model <-
function(name,name_index,nb_cores=1,seed=3,nb_tree=5000,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),raref=FALSE,empirical=TRUE,randomize=TRUE,nb_random=10,provided_tree=NULL,tolerance=0.05,overwrite=TRUE,figure=FALSE,path=NULL,path_alignment=NULL,...){
  
  if(!exists("name")) stop(print("Please provide the name of the dataset "))
  if(!exists("name_index")) stop(print("Please provide the name of the different OTU alignments "))
  if(is.null(path)) {path <- getwd()}
  if(!is.character(path)) {path <- getwd()}
  setwd(path)
  
  if (is.null(path_alignment)){ path_alignment <- path}
  
  if (is.null(provided_tree)){
    if (!file.exists(paste("host_tree_",name,".tre",sep=""))){stop(print("Please provide the host tree (format .tre) in the working directory"))}
    provided_tree <- read.tree(paste("host_tree_",name,".tre",sep=""))
  }
  host_tree <- provided_tree
  if (!is.binary(host_tree)) stop(print("Please provide a binary host tree"))
  if (!is.rooted(host_tree)) stop(print("Please provide a rooted host tree"))
  if (!is.ultrametric(host_tree)) stop(print("Please provide an ultrametric host tree"))
  
  print("Data preparation:")
  output <- mclapply(1:length(name_index), prepare_data_HOME, mc.cores=nb_cores,name=name,name_index=name_index,provided_tree=host_tree,path=path,path_alignment=path_alignment)
  
  print("Tree bank simulations:")
  output <- mclapply(1:length(lambda),simul_bank_tree,mc.cores=nb_cores,name=name,provided_tree=host_tree,nb_tree=nb_tree,lambda=lambda,seed=seed)
  
  print("Global inference:")
  for (index in name_index){output <- fit_HOME(index=index,name=name,nb_tree=nb_tree,lambda=lambda,nb_cores=nb_cores,raref=raref,tolerance=tolerance)}
  
  print("Initial output:")
  output <- mclapply(1:length(name_index), output_results_HOME, mc.cores=nb_cores,name=name,name_index=name_index,lambda=lambda,nb_tree=nb_tree,empirical=empirical,randomize=F,raref=raref, figure=figure)
  
  if (randomize==T){
    print("Model selection:")
    for (index in name_index){output <- model_selection_HOME(index=index,name=name,nb_tree=nb_tree,lambda=lambda,nb_cores=nb_cores,seed=seed,nb_random=nb_random,overwrite=overwrite,tolerance=tolerance)}
    
    print("Output:")
    output <- mclapply(1:length(name_index), output_results_HOME, mc.cores=nb_cores,name=name,name_index=name_index,lambda=lambda,nb_tree=nb_tree,empirical=empirical,randomize=T,raref=raref, figure=figure)
  }
}
