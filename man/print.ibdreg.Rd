% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:22:58 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.ibdreg.Rd,v 1.3 2010/02/08 18:22:58 sinnwell Exp $ -->
% $Locker:  $ -->
% 
% Log: print.ibdreg.sgml,v $
% Revision 1.5  2008/10/14 20:46:51  sinnwell
% addkeywords
%
\name{print.ibdreg}
\alias{print.ibdreg}
\title{
  print and ibdreg object
}
\description{
Print results stored in an ibdreg object, showing the function call,
the counts of relative pairs used, and the tests performed for linkage
with and without covariates on the relative pairs. A table is given to
show the test at the chromosome position which had the minimum p-value.
}
\usage{
\method{print}{ibdreg}(x, digits=max(options()$digits - 2, 5), ...)
}
\arguments{

\item{x}{
An ibdreg object
}
\item{digits}{
number of significant digits to print for numeric values
}
\item{\dots }{
optional parameters for the print method
}
}
\value{
nothing is returned
}
\seealso{
\code{\link{ibdreg}}
}
% docclass is function
% Converted by Sd2Rd version 43267.
