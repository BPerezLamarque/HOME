\name{simulate_data}
\alias{simulate_data}

\title{
Simulate microbial alignments according to different models of Host-Microbiota evolution
}

\description{
Runs the model inference for all the different symbiont (called index) on the same host tree.

}

\usage{
simulate_data(name, mu=1, n=20, seed=1,nb_cores=1,name_index,
simul,N=300,proportion_variant=0.1,path,provided_tree,...)}

\arguments{
\item{name}{
the name of the run
}

\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs) to simulate
}
\item{path}{
the path to the working directory
}

\item{nb_cores}{
a number of cores to run the analyses
}

\item{seed}{
a seed to assure the reproducibility
}

\item{simul}{
a vector of scenarii for each simulated alignment ("0" for strict vertical transmission, any positive integer for the number of host-switches transmissions, and "indep" to simulate independent evolutions and environmental acquisitions). This vector must have the same length as "name_index"
}
\item{N}{
a positive interger value representing the total number of nucleotides (variant and invariant nucleotides)
}
\item{proportion_variant}{
a value between 0 and 1 representing the proportion of variant nucleotides among N
}

\item{provided_tree}{
optional - the n-tips, binary, rooted, and ultrametric host tree to simulate the alignments on it (by default a pure-birth host tree is simulated)}

}


\value{
The function simulate microbial nucleotidic alignments according to different models of Host-Microbiota evolution (strict vertical transmission on the host tree, )

}
\details{

If you provide a binary rooted ultrametric host tree, it must be in Newick format (.tre).
}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Modeling host-microbiota evolution for characterizing symbiont inheritance – application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{HOME_model}}
}
\examples{

# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)
simulate_data(name="example_simulation", mu=1, n=10,name_index=c("Simul_1","Simul_2","Simul_3"), simul=c(0,5,"indep"), N=300,
proportion_variant=0.1)


}





























