\name{resources}
\alias{resources}
\title{
  track system resources
}
\description{
track system resources during evaluation of an expression
}
\usage{
resources(expr, doPrint=TRUE)
}
\arguments{
\item{expr}{
An expression to be evaluated by eval() within the sys.parent()
evaluation frame.
}
\item{doPrint}{
Logical, print the resource information in the call?
}
}
\value{
  Return the data frame that is printed by default if the object is
  saved to a value (via invisible()). The expression is evaluated and
  details are printed for CPU time, overall time, child process time,
  and heap space memory usage.
}
\details{
The resources given by Venebles and Ripley had reported both Cache and
Working memory usage.  S-PLUS 7 no longer tracks Cache usage, so we
only report the Working.  Similarly, R has Ncells and Vcells, and the
more important measure for users is Vcells, the memory heap.
}
\section{References}{
Venebles WN, Ripley BD. "Statistics and Computing."
Springer-Verlag, New York, NY, 2000.  (page 151.)
}
\seealso{
\code{\link{proc.time}},
R: gc, S-PLUS: mem.tally.reset, mem.tally.report
}
\examples{
resources({
  norm.dat = rnorm(10000) 
  norm.mat = matrix(rnorm(10000), nrow=100)
})
}
% docclass is function
% Converted by Sd2Rd version 43267.
