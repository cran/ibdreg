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
% Converted by Sd2Rd version 1.21.
