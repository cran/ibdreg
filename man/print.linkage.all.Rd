% $Author: sinnwell $ -->
% $Date: 2008/10/15 16:30:32 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.linkage.all.Rd,v 1.2 2008/10/15 16:30:32 sinnwell Exp $ -->
% $Locker:  $ -->
%
%$Log: print.linkage.all.Rd,v $
%Revision 1.2  2008/10/15 16:30:32  sinnwell
%usage
%
%Revision 1.1  2008/10/14 21:52:10  sinnwell
%Initial revision
%
%Revision 1.5  2008/10/14 20:47:07  sinnwell
%addkeywords
%
\name{print.linkage.all}
\alias{print.linkage.all}
\title{
  Print or summarize a linkage.all object
}
\description{
Print a linkage.all object, made within ibdreg. For the linkage test,
print the test statistic, degrees of freedom, and pvalue at the
position with the smallest pvalue.
}
\usage{
\method{print}{linkage.all}(x, digits=max(options()$digits - 2, 5), ...)
}
\arguments{
  
\item{x}{
A linkage.all object
}
\item{digits}{
number of significant digits to print for numeric values
}
\item{\dots }{
additional display parameters
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
\code{\link{linkage.all}},
\code{\link{plot.linkage.all}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
