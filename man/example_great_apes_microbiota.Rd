\name{HOME_model}
\alias{HOME_model}

\title{
Example of empirical application of HOME
}

\description{
This function provide an example of empirical application using HOME for 3 different OTUs from the great apes microbiota. It only download the data (alignments and host tree).
}

\usage{example_great_apes_microbiota(name,path)}

\arguments{
\item{name}{
the name of the run
}

\item{path}{
the path to the working directory
}
}

\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Modeling host-microbiota evolution for characterizing symbiont inheritance – application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{simulate_data}}
}
\examples{
# This example may take a little bit of time. Be patient!

example_great_apes_microbiota(name="great_apes",path=getwd())


# Run HOME (parameters estimations and model selections)

HOME_model(name="great_apes",name_index=c("OTU892624276","OTU47610657","OTU733943228"), nb_tree=500,
lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25), empirical=TRUE,randomize=TRUE,nb_random=10,path)

# NB: this should be run on a larger number of tree!

}




























