\name{prepare_data}
\alias{prepare_data}

\title{
Prepare the data for the inference
}

\description{
Prepare the data for the inference (format, model of nucleotidic evolution...).

}

\usage{prepare_data(iter,name,name_index,path,...)
}

\arguments{
\item{name}{
the name of the run
}
\item{iter}{
the number of iteration of a particular index from name_index
}
\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs)
}
\item{path}{
the path to the working directory
}

}


\value{
The function creates files to run the following functions simul_bank_tree, global_inference...
}
\details{
Please provide the binary rooted ultrametric host tree (Newick format .tre) in the working directory. Its name must be weel-formated "host_tree_NAME.tre" (match the name of the run). A folder (“path_alignment” or "path" by default") must contain all the alignments with the filenames "alignment_NAME_OTU.fas" (FASTA format).
See reference for more details.

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Modeling host-microbiota evolution for characterizing symbiont inheritance – application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{simulate_data}}
}
\examples{
# Some examples may take a little bit of time. Be patient!


# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)

name="example_simulation"
name_index=c("Simul_1","Simul_2","Simul_3")
path=getwd()

simulate_data(name, mu=1, n=10,name_index, simul=c(0,5,"indep"), N=300, proportion_variant=0.1)


# Inference

# Prepare the data (format, substitution model...)
for (i in 1:3){prepare_data(iter=i,name,name_index,path)}

# Simulate the bank of trees
for (ksi in seq(1,15)){simul_bank_tree(ksi,name,nb_tree=1000,lambda=seq(1,15),seed=1)}

# Infer the parameters
for (i in 1:3){global_inference(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1)}

# Plot the first outputs
for (i in 1:3){output_results(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
for (i in 1:3){model_selection(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1,seed=1)}

# Plot the final outputs
for (i in 1:3){output_results(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=TRUE,raref=FALSE)}



}





























