% $Author: sinnwell $ -->
% $Date: 2008/10/14 19:39:15 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/exact.ibd.var.Rd,v 1.1 2008/10/14 19:39:15 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: exact.ibd.var.Rd,v $
%Revision 1.1  2008/10/14 19:39:15  sinnwell
%Initial revision
%
%Revision 1.2  2006/11/21 17:23:08  sinnwell
%usage was still create.ibd.var, change to exact.ibd.var
%
%Revision 1.1  2006/10/12 14:03:13  sinnwell
%Initial revision
%
%Revision 1.3  2006/10/10 21:19:38  sinnwell
%format return object
%
%Revision 1.2  2006/08/07 17:44:03  sinnwell
%add commented example, add more details about perl script
%
%Revision 1.1  2006/04/19 22:15:00  sinnwell
%Initial revision
%
\name{exact.ibd.var}
\alias{exact.ibd.var}
\title{
  create an ibd.var object
}
\description{
create an ibd.var object from a temporary output file produced by
exact.ibd.var.pl perl script, which uses exact computations from merlin 
}
\usage{
exact.ibd.var(file)
}
\arguments{
\item{file}{
a temporary output file created by exact.ibd.var.pl, a perl scripts
provided in /ibdreg/perl/ in the ibdreg package
}

}
\value{
an ibd.var object (e.g. ret) that contains the following elements for each
pedigree:
<pre>
ret\$ped.id: pedigree id

ret\$person1.id: vector of ids for first person in the relative pair

ret\$person2.id: vector of ids for second person in the relative pair

ret\$sm: mean vector of ibd sharing between relative pairs \{person1.id, person2.id\}

ret\$sv: variance-covariance matrix for ibd sharing between pairs of
relative pairs.
</pre
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{sim.ibd.var}}
}
\examples{
## create a temporary file using perl script

# unix\% exact.ibd.var.pl chrom1.pre 1 chr1.var.tmp

## make an ibd.var object from chr1.var.tmp file

# RorS> chr1.ibd.var <- exact.ibd.var("chr1.var.tmp")
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
