% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:19:59 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.ibdreg.unexpect.Rd,v 1.3 2010/02/08 18:19:59 sinnwell Exp $ -->
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
\seealso{
\code{\link{plot.ibdreg}},
\code{\link{ibdreg}}
}
% docclass is function
% Converted by Sd2Rd version 43267.
