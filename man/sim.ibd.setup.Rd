% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:39:33 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.ibd.setup.Rd,v 1.2 2010/02/08 18:39:33 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: sim.ibd.setup.Rd,v $
%Revision 1.2  2010/02/08 18:39:33  sinnwell
%rm empty sections
%
%Revision 1.1  2008/10/14 18:09:23  sinnwell
%Initial revision
%
%Revision 1.2  2006/10/27 17:20:18  sinnwell
%fix ped description
%
%Revision 1.1  2006/10/10 21:20:38  sinnwell
%Initial revision
%
\name{sim.ibd.setup}
\alias{sim.ibd.setup}
\alias{setFounderAlleles}
\title{
  set-up functions for sim.mark.prop
}
\description{
set-up functions for sim.mark.prop
}
\usage{
sim.ibd.setup(ped, miss.val = c(NA, 0), x.linked=FALSE)
setFounderAlleles(ped, x.linked=FALSE)
}
\arguments{

\item{ped}{
a data.frame with data to characterize one pedigree.  The columns
contain id codes for person, mother, and father; also a column for sex.
}
\item{miss.val}{
code(s) for missing values ped
}
\item{x.linked }{
logical; if TRUE, apply method to x chromosome
}
}
\value{
A ped object to be used in sim.mark.prop
}
\seealso{
\code{\link{sim.mark.prop}}
}
% docclass is function
% Converted by Sd2Rd version 43267.
