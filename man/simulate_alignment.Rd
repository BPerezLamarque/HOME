\name{simulate_alignment}
\alias{simulate_alignment}

\title{
Simulate a sequence alignment on a phylogentic tree
}

\description{
This function simulates a sequence alignment on a given phylogenetic tree according to a K80 model of DNA substitution.
}

\usage{
simulate_alignment(tree, mu=0.5, N=300, proportion_variant=0.1)
}

\arguments{

\item{tree}{
a n-tips, binary, rooted, and ultrametric tree to simulate the alignments on it. All branch lengths must be >0. }


\item{mu}{
a positive value representing the simulated substitution rate.
}

\item{N}{
a positive integer value representing the total number of nucleotides (variable and unvariable nucleotides)
}

\item{proportion_variant}{
a value between 0 and 1 representing the proportion of variable nucleotides among N (i.e. each site has a probability "proportion_variant" to be variable)
}

}


\value{
The function outputs a nucleotidic alignment (as a matrix with each row corresponding to a sequence per tip of the tree).

}

\references{
Perez-Lamarque B, Morlon H (2019). Characterizing symbiont inheritance during host-microbiota evolution: Application to the great apes gut microbiota. Molecular Ecology Resources 19:1659-1671.
}
\author{Benoit Perez-Lamarque}
\seealso{
\code{\link{sim_microbiota}}
\code{\link{HOME_model}}
}
\examples{

set.seed(3)
tree <- phytools::pbtree(n=20) # simulate host phylogeny


# simulate an alignment
simulate_alignment(tree, mu=0.1, N=300, proportion_variant=0.1)


}


