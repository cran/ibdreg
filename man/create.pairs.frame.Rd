\name{create.pairs.frame}
\alias{create.pairs.frame}
\alias{create.pairs.frame.cvec}
\title{
  create a data frame for relative pairs
}
\description{
create a data frame for relative pairs containing identifiers,
response, and covariates
}
\usage{
create.pairs.frame(ibd.dat, model.mat, formula, c.scale)
create.pairs.frame.cvec(c.scale, prior2, prior1)
}
\arguments{

\item{ibd.dat }{
object of class ibd.dat
}
\item{model.mat }{
a model matrix containing covariates for each person; to be used in
making covariates for the relative pair
}
\item{formula}{
format is either ~1 or ~pairs.fun(cov1).  If the latter, the columns 
of the returned x.mat will be the result of the function operated on
the covariate for each relative in the relative pair.
}
\item{c.scale }{
scaling factor for covariates.  Either "nodom" for no dominance 
variance or "minimax".
}
\item{prior2}{
Prior probability of relative pairs to share 2 alleles ibd (from
ibd.dat object)
}
\item{prior1}{
Prior probability of relative pairs to share 1 allele ibd (from
ibd.dat object)
}
}
\value{
a list of 3 data frames all with the same number of rows.  These data
frames are id.df (3-columns of ped.id, person1.id, person2.id), y.mat
(estimated ibd allele sharing, which is the regression response), status.df
(3-columns for AA, UU, or AU status for the pair), x.mat (model matrix
with intercept and, if covariates in formula, columns for covariates).
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{create.ibd.dat}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
