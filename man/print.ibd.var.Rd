% $Author: sinnwell $ -->
% $Date: 2008/10/15 16:30:32 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.ibd.var.Rd,v 1.2 2008/10/15 16:30:32 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: print.ibd.var.Rd,v $
%Revision 1.2  2008/10/15 16:30:32  sinnwell
%usage
%
%Revision 1.1  2008/10/14 20:49:41  sinnwell
%Initial revision
%
%Revision 1.1  2006/11/07 20:04:15  sinnwell
%Initial revision
%

\name{print.ibd.var}
\alias{print.ibd.var}
\title{
  print an ibd.var object
}
\description{
Print results stored in an ibd.var object.  
If ped.id is NULL and sinkfile is NULL, print everything to a made-up
sink file because the output is most likely too big for standard
output.  If ped.id is not NULL, print the ibd.var elements that 
match ped.ids.  If sinkfile is not NULL, sink results to that file.
}
\usage{
 \method{print}{ibd.var}(x, ped.id=NULL, sinkfile=NULL, digits=max(options()$digits - 2, 5), ...)
}
\arguments{

\item{x}{
An ibd.var object, created by either exact.ibd.var or sim.ibd.var
}
\item{ped.id }{
Vector of pedigree IDs to specify a subset of the pedigrees to print
}
\item{sinkfile}{
Name of file to sink too.  Printed output from ibd.var is large,
especially if ped.id is not given, so enforce a file name for where to
sink print results.  If ped.id is NULL and sinkfile is not given,
default is set to name of x with ".sink".
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
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{sim.ibd.var}},
\code{\link{exact.ibd.var}},
\code{\link{ibdreg}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
