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
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{sim.mark.prop}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
