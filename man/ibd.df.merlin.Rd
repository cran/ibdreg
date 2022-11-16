% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:10:06 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/ibd.df.merlin.Rd,v 1.2 2010/02/08 18:10:06 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: ibd.df.merlin.Rd,v $
%Revision 1.2  2010/02/08 18:10:06  sinnwell
%rm empty section
%
%Revision 1.1  2008/10/14 19:39:03  sinnwell
%Initial revision
%
%Revision 1.5  2008/10/14 18:37:34  sinnwell
%keywords
%
\name{ibd.df.merlin}
\alias{ibd.df.merlin}
\title{
  Create and return a data.frame with ibd data from merlin
}
\description{
Create and return a dataframe of ped.id, person1.id, 
person2.id, and ibd information from merlin
}
\usage{
ibd.df.merlin(ibd.dat)
}
\arguments{

\item{ibd.dat }{
data.frame with posterior IBD data from Merlin, with column 
names: FAMILY, ID1, ID2, MARKER
}
}
\value{
The returned data.frame, call it df, with the following elements:
  \item{ped.id}{pedigree id}
  \item{person1.id}{person 1 of a relative pair}
  \item{person2.id}{ person 2 identifier of a relative pair within a pedigree}
  \item{post0}{matrix objects that can be referenced from df.  Each}
  \item{post1}{is a matrix with posterior ibd sharing probabilities for}
  \item{post2}{relative pairs (rows) at each chromosome position (columns)}
}
