\name{mergeIBD}
\alias{mergeIBD}
\title{
  Merge ibd and covariate data
}
\description{
Merge ibd probabilities data frame and an individual's sex 
on pedigree id and person id
}
\usage{
mergeIBD(ibd.dat, sex.dat)
}
\arguments{

\item{ibd.dat }{
A data frame with ped.id, person1.id, person2.id, and the
estimated ibd sharing values at each chromosome position for the two people.
}
\item{sex.dat }{
a data frame containing ped.id, person.id, and sex
}
}
\value{
a data frame with identifiers for relative pairs, their ibd data, and
the two sex codes for the pair
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{

}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
