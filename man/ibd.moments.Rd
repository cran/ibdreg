% $Author: sinnwell $ 
% $Date: 2010/02/08 18:10:15 $ 
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/ibd.moments.Rd,v 1.2 2010/02/08 18:10:15 sinnwell Exp $ -->
% $Locker:  $
% 
%$Log: ibd.moments.Rd,v $
%Revision 1.2  2010/02/08 18:10:15  sinnwell
%rm empty sections
%
%Revision 1.1  2008/10/14 18:05:13  sinnwell
%Initial revision
%
%Revision 1.1  2006/10/10 21:20:05  sinnwell
%Initial revision
%

\name{ibd.moments}
\alias{ibd.moments}
\title{
  calculate moments (mean, variance) for simulated ibd vectors
}
\description{
calculate moments (mean, variance) for simulated ibd vectors
for an individual predigree from sim.mark.prop
}
\usage{
ibd.moments(simped, person.indx, male, x.linked=FALSE)
}
\arguments{

\item{simped}{
The object returned from sim.mark.prop
}
\item{person.indx }{
person identifier
}
\item{male}{
binary code (1/0) for if person is male
}
\item{x.linked }{
logical, if TRUE, apply the method for an x-linked chromosome
}
}
\value{
a list with person1.id person2.id sm and sv
}
\seealso{
\code{\link{sim.mark.prop}}
}
% docclass is function
% Converted by Sd2Rd version 43267.
