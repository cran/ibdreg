% $Author: sinnwell $ -->
% $Date: 2008/10/15 16:28:17 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.linkage.tests.Rd,v 1.2 2008/10/15 16:28:17 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: plot.linkage.tests.Rd,v $
%Revision 1.2  2008/10/15 16:28:17  sinnwell
%usage
%
%Revision 1.1  2008/10/14 20:43:43  sinnwell
%Initial revision
%
%Revision 1.3  2008/10/14 20:30:52  sinnwell
%addkeywords
%

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
\method{plot}{linkage.tests}(x, ...)
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
% Converted by Sd2Rd version 43267.
