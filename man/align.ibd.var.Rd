% $Author: sinnwell $ -->
% $Date: 2008/10/14 19:25:10 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/align.ibd.var.Rd,v 1.1 2008/10/14 19:25:10 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: align.ibd.var.Rd,v $
%Revision 1.1  2008/10/14 19:25:10  sinnwell
%Initial revision
%
%Revision 1.3  2008/10/14 18:36:46  sinnwell
%add keywords
%

\name{align.ibd.var}
\alias{align.ibd.var}
\title{
  align ibd.var object with relative pairs within pairs.frame object
}
\description{
align ibd.var object with relative pairs within pairs.frame object
}
\usage{
align.ibd.var(id.df, ibd.var, epsilon=1e-5)
}
\arguments{
\item{id.df }{
a data.frame containing the three id columns (ped.id, person1.id,
person2.id) from a pairs.frame object
}
\item{ibd.var }{
an ibd.var object
}
\item{epsilon}{
cutoff for singular values in a generalized inverse (Ginv)
}
}
\value{
a list resembling an ibd.var object
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{

}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
