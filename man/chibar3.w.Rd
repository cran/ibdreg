% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:03:57 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/chibar3.w.Rd,v 1.2 2010/02/08 18:03:57 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: chibar3.w.Rd,v $
%Revision 1.2  2010/02/08 18:03:57  sinnwell
%rm empty sections
%
%Revision 1.1  2008/10/14 19:39:42  sinnwell
%Initial revision
%
%Revision 1.3  2008/10/14 18:36:55  sinnwell
%add keywords
%
\name{chibar3.w}
\alias{chibar3.w}
\title{
  Compute mixture proportions, w, for mixture chi-square
distribution
}
\description{
Compute mixture proportions, w, for mixture chi-square
distribution when k=0..3, where k = df of chi-squares
}
\usage{
chibar3.w(sigma)
}
\arguments{

\item{sigma}{
covariance matrix (3 x 3) of the 3 distributions
}
}
\value{
mixing proportions, with element w[k+1] the mixture proportion 
for k, k=0..3.
}
\seealso{
\code{\link{pchibar}}
}
\examples{
chibar3.w(sigma=diag(c(1,1,1)))
}
% docclass is function
% Converted by Sd2Rd version 43267.
