% $Author: sinnwell $ -->
% $Date: 2008/10/15 16:30:32 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.linkage.tests.Rd,v 1.2 2008/10/15 16:30:32 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: print.linkage.tests.Rd,v $
%Revision 1.2  2008/10/15 16:30:32  sinnwell
%usage
%
%Revision 1.1  2008/10/14 21:58:59  sinnwell
%Initial revision
%
%Revision 1.4  2008/10/14 20:47:31  sinnwell
%addkeywords
%

\name{print.linkage.tests}
\alias{print.linkage.tests}
\title{
  Print a linkage.tests object
}
\description{
Print a linkage.tests object, made within ibdreg. For each test, print
the score test statistic, degree(s) of freedom, and pvalue at the
position with the smallest pvalue.
}
\usage{
\method{print}{linkage.tests}(x, digits=max(options()$digits - 2, 5), 
                    show.model.tests=FALSE, ...)
}
\arguments{
\item{x}{
A linkage.tests object
}
\item{digits}{
The number of significant digits to print for numeric values
}
\item{show.model.tests }{
logical, if TRUE, show tests with model-based covariance, in addition
to the other tests
}
\item{\dots }{
Additional print parameters
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
\code{\link{linkage.tests}},
\code{\link{plot.linkage.tests}}
}
\examples{
# reg.out <- ibdreg(formula=~1, ped.id=ped.id, person.id=person.id, 
#                   status=status, status.method="AA", data=data,
#                   ibd.dat=ibd.dat.obj, ibd.var=ibd.var.obj) 

# print.linkage.tests(reg.out$AA.linkage, digits=4)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
