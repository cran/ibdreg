% $Author: sinnwell $ -->
% $Date: 2008/10/15 16:27:21 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.ibdreg.Rd,v 1.2 2008/10/15 16:27:21 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: $
%
\name{plot.ibdreg}
\alias{plot.ibdreg}
\title{
  plot an ibdreg object
}
\description{
plot test results within an ibdreg object
}
\usage{
 \method{plot}{ibdreg}(x, ...)
}
\arguments{

\item{x}{
an ibdreg object
}
\item{\dots }{
optional parameters for plot function, title is expected
}
}
\value{
nothing is returned
}
\section{Side Effects}{

}
\details{
If status.method was AA, UU, or AU, only tests for linkage on
the chosen group of relative pairs will be plotted.  If status.method was ALL,
then the linkage test in linkage.all object is plotted.  If
status.method was ALL and only linkage tests were done, 4 lines for
linkage of AA, UU, AU, ALL linkage are plotted together.
}
\section{References}{

}
\seealso{
\code{\link{ibdreg}},
\code{\link{plot.ibdreg.unexpect}},
\code{\link{plot.linkage.tests}},
\code{\link{plot.linkage.all}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
