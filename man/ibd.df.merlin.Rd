\name{ibd.df.merlin}
\alias{ibd.df.merlin}
\title{
  Create and return a data.frame with ibd data from merlin
}
\description{
Create and return a dataframe of ped.id, person1.id, 
person2.id, and ibd information from merlin
}
\usage{
ibd.df.merlin(ibd.dat)
}
\arguments{

\item{ibd.dat }{
data.frame with posterior IBD data from Merlin, with column 
names: FAMILY, ID1, ID2, MARKER
}
}
\value{
The returned data.frame, call it df, with the following elements:

df\$ped.id:     pedigree id

df\$person1.id:
df\$person2.id: identifiers for the relative pair within a pedigree

df\$post0:      matrix objects that can be referenced from df.  Each
df\$post1:      is a matrix with posterior ibd sharing probabilities for
df\$post2:      relative pairs (rows) at each chromosome position (columns)
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
