\name{chibar3.w}
\alias{chibar3.w}
\title{
  Compute mixture proportions, w, for mixture chi-square
distribution
}
\description{
Compute mixture proportions, w, for mixture chi-square
distribution when k=0..3, where k = df of chi-squares
}
\usage{
chibar3.w(sigma)
}
\arguments{

\item{sigma}{
covariance matrix (3 x 3) of the 3 distributions
}
}
\value{
mixing proportions, with element w[k+1] the mixture proportion 
for k, k=0..3.
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{pchibar}}
}
\examples{
chibar3.w(sigma=diag(c(1,1,1)))
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
