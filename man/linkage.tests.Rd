% $Author: sinnwell $ -->
% $Date: 2008/10/14 19:46:24 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.tests.Rd,v 1.1 2008/10/14 19:46:24 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: $
%
\name{linkage.tests}
\alias{linkage.tests}
\alias{linkage.tests.B0}
\alias{linkage.tests.pedloop}
\title{
  Linkage tests for relative pair groups
}
\description{
Linkage tests, with and without covariates, for relative pair status groups
Affected-Affected (AA), Unaffected-Unaffected (UU), or Affected-Unaffected (AU)
}
\usage{
linkage.tests(y.mat, x.adj, ibdvar.lst, status.method, epsilon=1e-5)
linkage.tests.B0(y.mat, xvec, ibdvar.lst)
linkage.tests.pedloop(y.vec, x.adj, ibdvar.lst, 
                      B0, concordant=TRUE, epsilon)
}
\arguments{

\item{y.mat }{
matrix of ibd sharing statistics, which are the response in the
regression model; each column represents one chromosome position.
}
\item{y.vec }{
vector of ibd sharing statistics (at one position), which are the response in the
regression model (linkage.tests.pedloop only)
}
\item{x.adj }{
matrix of covariates, x, adjusted by c.scale.
}
\item{xvec}{
the intercept column from x.adj, needed for pre-calculating Beta-0
}
\item{ibdvar.lst }{
a list containing an element for each pedigree. The elements include
ped.id, person1.id, person2.id, and the variance-covariance matrix 
of pairs of pairs of subjects.
}
\item{status.method }{
Character string indicating which relative pairs to apply the linkage
tests: "AA", "UU", or "AU".
}
\item{B0}{
The intercept in the regression model, ("beta-zero").
}
\item{concordant}{
logical; indicate whether to perform tests for linkage on
concordant pairs (AA, UU), or discordant (AU).
}
\item{epsilon}{
cutoff for singular values in generalized inverse (Ginv)
}
}
\value{
an object of class linkage.tests
}
\section{Side Effects}{

}
\details{
The calculations for linkage tests require a pre-calculation of
beta-zero (linkage.tests.B0) and then proceed with calculating the
linkage test statistics for each position by looping over all
pedigrees (linkage.tests.pedloop).  Linkage tests can only be applied
to one status group at a time, because the tests are restricted to
whether the relative pairs are all concordant or discordant.  
Unexpected allele sharing can be detected when unconstrained 
linkage > constrained linkage tests.  The method also returns z scores
for pedigree-specific contributions to unexpected allele sharing
indicated in linkage tests.
}
\section{References}{

}
\seealso{
\code{\link{ibdreg}},
\code{\link{print.linkage.tests}},
\code{\link{plot.linkage.tests}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
