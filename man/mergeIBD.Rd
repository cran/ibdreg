% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:14:11 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/mergeIBD.Rd,v 1.2 2010/02/08 18:14:11 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: $
%
\name{mergeIBD}
\alias{mergeIBD}
\title{
  Merge ibd and covariate data
}
\description{
Merge ibd probabilities data frame and an individual's sex 
on pedigree id and person id
}
\usage{
mergeIBD(ibd.dat, sex.dat)
}
\arguments{

\item{ibd.dat }{
A data frame with ped.id, person1.id, person2.id, and the
estimated ibd sharing values at each chromosome position for the two people.
}
\item{sex.dat }{
a data frame containing ped.id, person.id, and sex
}
}
\value{
a data frame with identifiers for relative pairs, their ibd data, and
the two sex codes for the pair
}
% docclass is function
% Converted by Sd2Rd version 43267.
