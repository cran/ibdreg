% $Author: sinnwell $ -->
% $Date: 2008/10/14 20:07:38 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/minpRows.Rd,v 1.1 2008/10/14 20:07:38 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: $
%
\name{minpRows}
\alias{minpRows}
\title{
  find the rows of a data frame with the minimum pvalue
}
\description{
find the rows of a data frame with the minimum pvalue, best used in a
print or summary function to summarize a group of test statistics.
}
\usage{
minpRows(obj, colnames=NULL, rowname=NULL, col.indx=ncol(obj))
}
\arguments{

\item{obj}{
a data frame for a group of statistical tests, with pvalues
}
\item{colnames}{
new names to give to columns
}
\item{rowname}{
new names to give the row(s).  If multiple rows contain the minimum
pvalue, rowname is assigned to each of them, appending a digit at the
end (1, 2, etc).
}
\item{col.indx }{
the index of the columns that contains the pvalue, default is the last
column of obj
}
}
\value{
a subset of the obj data.frame rows
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{

}
\examples{
## create a data frame of chi-square(1) tests
tests <- data.frame(pos=1:3, chitest=c(3,4,2), df=rep(1,3), 
                    pval=1-pchisq(c(3,4,2), df=1))

# find the rows of the best tests
best.tests <- minpRows(tests, colnames=c("position", "chi.test", "df","pvalue"))
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 43267.
