\name{plot.linkage.tests}
\alias{plot.linkage.tests}
\title{
  Plot a linkage.tests object
}
\description{
Plot a linkage tests object, which was made within ibdreg.  Plot the
-log10(pvalue) for pvalues from the score test statistics for linkage
with and without covariates, constrained for one-sided tests.
}
\usage{
plot.linkage.tests(x, ...)
}
\arguments{

\item{x}{
An object of class linkage.tests
}
\item{\dots }{
additional plot parameters, this method checks for 'title'
}
}
\value{
nothing is returned
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{ibdreg}},
\code{\link{linkage.tests}}
}
\examples{
# reg.out <- ibdreg(formula=~1, ped.id=ped.id, person.id=person.id, 
#                   status=status, status.method="AA", data=data,
#                   ibd.dat=ibd.dat.obj, ibd.var=ibd.var.obj) 

# plot.linkage.tests(reg.out$AA.linkage, title="AA, formula: ~1")
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
