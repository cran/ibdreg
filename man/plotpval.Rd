% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:20:50 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plotpval.Rd,v 1.2 2010/02/08 18:20:50 sinnwell Exp $ -->
% $Locker:  $ -->
%
%$Log: plotpval.Rd,v $
%Revision 1.2  2010/02/08 18:20:50  sinnwell
%rm empty section
%
%Revision 1.1  2008/10/15 14:03:57  sinnwell
%Initial revision
%
%Revision 1.5  2008/10/14 20:46:33  sinnwell
%addkeywords
%

\name{plotpval}
\alias{plotpval}
\title{
  Plot the -log10 of pvalues
}
\description{
Plot lines connecting the the -log10 of pvalues for a series of tests
}
\usage{
plotpval(pos, pmat, lty=1, lwd=1, col=1, title="", legend, 
         xlab="Position(cM)", ylab="-log10(pvalue)")
}
\arguments{
\item{pos}{
an index for the series of tests, usually assuming adjacent tests are
correlated.
}
\item{pmat}{
matrix of p-values, each column represents a different test, rows
correspond to the pos[ition] index
}
\item{lty}{
line type, see par()
}
\item{lwd}{
line width, see par()
}
\item{col}{
line color, see par()
}
\item{title}{
title for the top-center of the plot
}
\item{legend}{
vector of text strings explaning the tests in each column of pmat
}
\item{xlab}{
x axis label
}
\item{ylab}{
y axis label
}
}
\value{
nothing is returned
}
\section{Side Effects}{
none
}
\details{
The default height for the y-axis is 3, which corresponds to a p-value
of 1e-3.  If any p-value is smaller, the height is -log10 of the
smallest p-value in pmat.
}
\seealso{
\code{\link{par}},
\code{\link{title}},
\code{\link{legend}}
}
\examples{
plotpval(pos=1:10, pmat=cbind(runif(10)/5, runif(10)/10), 
         lty=c(1,2), col=c(1,2), legend=c("runif/5", "runif/10"))
}
% docclass is function
% Converted by Sd2Rd version 43267.
