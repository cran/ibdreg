% $Author: sinnwell $ -->
% $Date: 2008/10/14 20:07:25 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/pairSum.Rd,v 1.1 2008/10/14 20:07:25 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: pairSum.Rd,v $
%Revision 1.1  2008/10/14 20:07:25  sinnwell
%Initial revision
%
%Revision 1.2  2006/10/30 19:58:18  sinnwell
%trying to make usage show both functions
%
%Revision 1.1  2006/09/12 21:26:02  sinnwell
%Initial revision
%
%Revision 1.2  2006/09/12 20:41:45  sinnwell
%remove reference
%
%Revision 1.1  2006/04/26 15:22:48  sinnwell
%Initial revision
%

\name{pairSum}
\alias{pairSum}
\alias{pairDiff}
\title{
  Combine person-specific covariates to pair-specific
}
\description{
Functions to combine person-specific covariates to pair-specific
}
\usage{
pairSum(cov1, cov2)
pairDiff(cov1, cov2)
}
\arguments{

\item{cov1}{

}
\item{cov2}{

}
}
\value{
Either the sum or the difference (absolute) of the two covariates
}
\section{Side Effects}{

}
\details{
For use within the ibdreg package to make covariates that are a function
of the two covariates in a relative pair. These functions are to be
used in the formula parameter of ibdreg, and the pair-specific
covariates can be expressed as pairSum(cov) or pairSum(cov, cov). 
}
\section{References}{

}
\seealso{
\code{\link{ibdreg}}
}
\examples{
pairSum(20,30)

pairDiff(20,30)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
