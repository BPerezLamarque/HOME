\name{simul_bank_tree}
\alias{simul_bank_tree}

\title{
Simulate a bank of trees
}

\description{
This function simulates the bank of microbial trees, used for the inferences, for a given number of switches.

}

\usage{simul_bank_tree(nb_ksi,name,nb_tree,lambda,seed=1,...)}

\arguments{
\item{name}{
the name of the run
}
\item{nb_ksi}{
a positive integer which represent the number of host-switches to simulate on the host tree
}
\item{nb_tree}{
a number of tree to simulate
}
\item{lambda}{
a vector of integer values of number of switches to test during estimations
}
\item{seed}{
a seed to assure the reproducibility
}

}


\value{
The function outputs trees in the folder simulated_trees/ for each value of number of switches. Then, use the function \emph{global_inference} to infer the results.
}
\details{

The functions \emph{prepare_data} and \emph{simul_bank_tree} must be run before.
See reference for more details.

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Characterizing symbiont inheritance during host-microbiota evolution : application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{simulate_data}},
\code{\link{prepare_data}},
\code{\link{global_inference}},
\code{\link{output_results}},
\code{\link{model_selection}},
\code{\link{HOME_model}}
}
\examples{
# Some examples may take a little bit of time. Be patient!


# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)

name="example_simulation"
name_index=c("Simul_1","Simul_2","Simul_3")
path=getwd()

simulate_data(name, mu=1, n=10, name_index, path=path, simul=c(0,5,"indep"), N=300, proportion_variant=0.1)

# Inference

# Prepare the data (format, substitution model...)
for (i in 1:3){prepare_data(iter=i,name,name_index,path)}

# Simulate the bank of trees
for (nb_ksi in 1:length(seq(1,15))){simul_bank_tree(nb_ksi,name,nb_tree=1000,lambda=seq(1,15),seed=1)}

# Infer the parameters
for (i in 1:3){global_inference(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1)}

# Plot the first outputs
for (i in 1:3){output_results(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
for (i in 1:3){model_selection(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1,seed=1)}

# Plot the final outputs
for (i in 1:3){output_results(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=TRUE,raref=FALSE)}



}





























