% $Author: sinnwell $ -->
% $Date: 2008/10/14 19:46:24 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.all.Rd,v 1.1 2008/10/14 19:46:24 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: $
%
\name{linkage.all}
\alias{linkage.all}
\alias{linkage.all.pedloop}
\title{
  Linkage test statistics for ALL relative pairs
}
\description{
Calculate and store linkage test statistics for ALL relative pairs.
The tests include 
}
\usage{
linkage.all(y.mat, x.adj, ibdvar.lst, epsilon=1e-5)
linkage.all.pedloop(y.vec, x.adj, ibdvar.lst)
}
\arguments{
\item{y.mat }{
matrix of ibd sharing statistics, the response in the regression
model; each column represents one chromosome position.
}
\item{y.vec }{
vector of ibd sharing statistics (at one position), the response in the regression
model (pedloop only)
}
\item{x.adj }{
matrix of covariates, x, which are indicators for which status group
the relative pairs belong, adjusted by c.scale. 
}
\item{ibdvar.lst }{
a list containing an element for each pedigree. The elements include
ped.id, person1.id, person2.id, and the variance-covariance matrix 
of pairs of pairs osf subjects.
}
\item{epsilon}{
cutoff for singular values in generalized inverse (Ginv)
}
}
\value{
a linkage.all object, containing tests for linkage on all relative
pairs
}
\section{Side Effects}{

}
\details{
linkage.all sets up the linkage calculations for each position, and
the calculations are performed over all pedigrees in linkage.all.pedloop.
}
\section{References}{

}
\seealso{
\code{\link{ibdreg}},
\code{\link{print.linkage.all}},
\code{\link{plot.linkage.all}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
