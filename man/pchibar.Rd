% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:16:51 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/pchibar.Rd,v 1.2 2010/02/08 18:16:51 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: pchibar.Rd,v $
%Revision 1.2  2010/02/08 18:16:51  sinnwell
%rm empty sections
%
%Revision 1.1  2008/10/14 20:08:41  sinnwell
%Initial revision
%
%Revision 1.2  2008/10/14 19:49:36  sinnwell
%addkeywords
%

\name{pchibar}
\alias{pchibar}
\title{
  cumulative probability of a chibar distribution
}
\description{
cumulative probability of a chi-bar distribution (mixed chi-square)
}
\usage{
pchibar(x, df, wt)
}
\arguments{
\item{x}{
quantile of chi-bar distribution
}
\item{df}{
vector degrees of freedom
}
\item{wt}{
vector of mixing proportions (weights)
}

}
\value{
cumulative density of chi-bar distribution
}
\seealso{
\code{\link{pchisq}}
}
\examples{
pch = pchibar(3, df=c(0,1), wt=c(.5,.5))
1-pch  # compare to 1-pchisq(3,1)
}
% docclass is function
% Converted by Sd2Rd version 43267.
