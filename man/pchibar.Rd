\name{pchibar}
\alias{pchibar}
\title{
  cumulative probability of a chibar distribution
}
\description{
cumulative probability of a chi-bar distribution (mixed chi-square)
}
\usage{
pchibar(x, df, wt)
}
\arguments{
\item{x}{
quantile of chi-bar distribution
}
\item{df}{
vector degrees of freedom
}
\item{wt}{
vector of mixing proportions (weights)
}

}
\value{
cumulative density of chi-bar distribution
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{pchisq}}
}
\examples{
pch = pchibar(3, df=c(0,1), wt=c(.5,.5))
1-pch  # compare to 1-pchisq(3,1)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
