\name{model_selection}
\alias{model_selection}

\title{
Perform the model selections
}

\description{
This function performs the model selections of strict vertical transmission and independent evolution.

}

\usage{model_selection(index,name,nb_tree=1000,lambda,
nb_cores=1,seed=1,nb_random=10,...)
}

\arguments{
\item{name}{
the name of the run
}
\item{index}{
the name of the particular index from name_index (e.g. name of one OTU)
}
\item{nb_cores}{
a number of cores to run the analyses
}
\item{nb_tree}{
a number of tree for Monte-Carlo estimation of the number of switches (a low value will give inaccurate results whereas a high value will increase the computation time)
}
\item{lambda}{
a vector of integer values of number of switches to test during estimations
}
\item{seed}{
a seed to assure the reproducibility
}
\item{nb_random}{
a integer value of the number of randomizations (R) to perform during the model selection testing independent evolutions
}

}


\value{
The function outputs different files. Use the function \emph{output_results} to interpret them.}
\details{

The functions \emph{prepare_data}, \emph{simul_bank_tree}, and \emph{global_inference} must be run before.
See reference for more details.

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Characterizing symbiont inheritance during host-microbiota evolution : application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{simulate_data}},
\code{\link{prepare_data}},
\code{\link{simul_bank_tree}},
\code{\link{global_inference}},
\code{\link{output_results}},
\code{\link{HOME_model}}
}
\examples{
# Some examples may take a little bit of time. Be patient!


# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)

name="example_simulation"
name_index=c("Simul_1","Simul_2","Simul_3")
path=getwd()

simulate_data(name, mu=1, n=10,name_index, simul=c(0,5,"indep"), path=path, N=300, proportion_variant=0.1)

# Inference

# Prepare the data (format, substitution model...)
for (i in 1:3){prepare_data(iter=i,name,name_index,path)}

# Simulate the bank of trees
for (ksi in 1:length(seq(1,15))){simul_bank_tree(ksi,name,nb_tree=nb_tree,lambda=seq(1,15),seed=1)}

# Infer the parameters
for (i in 1:3){global_inference(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1)}

# Plot the first outputs
for (i in 1:3){output_results(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
for (i in 1:3){model_selection(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1,seed=1)}

# Plot the final outputs
for (i in 1:3){output_results(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=TRUE,raref=FALSE)}



}





























