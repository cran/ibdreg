% $Author: sinnwell $ -->
% $Date: 2008/10/15 16:27:21 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.ibdreg.unexpect.Rd,v 1.2 2008/10/15 16:27:21 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: $
%
\name{plot.ibdreg.unexpect}
\alias{plot.ibdreg.unexpect}
\title{
  Plot unconstrained and constrained linkage tests to
indicate areas of unexpected allele sharing
}
\description{
Plot the -log10(pvalue) of the unconstrained and constrained linkage
tests to indicate areas of unexpected allele sharing
}
\usage{
 \method{plot}{ibdreg.unexpect}(x, status.method="AA", ...)
}
\arguments{
\item{x}{
an ibdreg object
}
\item{status.method }{
Character string indicating which affection status for relative pairs
("AA", "AU", or "UU") is to be examined for unexpected IBD sharing.
Character string indicating which relative pairs to display linkage
tests for: "AA", "AU", or "UU".  These correspond to:
Affected-Affected, Affected-Unaffected, and Unaffected-Unaffected, respectively.
}
\item{\dots }{
extra parameters for the plot method.  This method looks for a
'title' parameter.
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
\code{\link{plot.ibdreg}},
\code{\link{ibdreg}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
