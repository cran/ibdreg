% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:12:35 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/ibdreg.Rd,v 1.2 2010/02/08 18:12:35 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%Log: ibdreg.sgml,v $
%revision 1.14  2008/10/14 18:37:47  sinnwell
%keywords
%
\name{ibdreg}
\alias{ibdreg}
\title{
  Regression Methods to Test for Linkage With Covariates
}
\description{
A method to test genetic linkage with covariates by regression models
that use the IBD status for relative pairs as the dependent variable
and pair-specific covariates as the independent variables.  Correlations of 
IBD statistics and covariates for relative pairs within the same
pedigree are accounted for.
}
\usage{
ibdreg(formula, status.method, c.scale='nodom', data,
                   status, ped.id, person.id,
                   ibd.dat, ibd.var, 
                   subset, weights, na.action,
                   min.pairs=1, epsilon=1e-5, ...)
}
\arguments{
\item{formula}{
either "~1" for no covariates to get only linkage tests, or a 
function of the paired covariates (e.g., "~pairs.sum(cov1, cov1)") for
tests that include covariates.
}
\item{status.method }{
Character string indicating which relative pairs to apply the
ibdreg method: "AA", "UU", "AU", or "ALL". 
These correspond to: Affected-Affected, Unaffected-Unaffected, 
Affected-Unaffected, and all pairs, respectively.
}
\item{ibd.dat }{
an ibd.dat object, created by create.ibd.dat() with elements: 
ped.id, person1.id, person2.id, which collectively identify relative
pairs, post0, post1, post2, which are the posterior probabilities of sharing 0,
1, and 2 alleles IBD at each chromosome position, and prior0, prior1, prior2, 
which are the null prior probabilities of sharing 0, 1, and 2 alleles ibd. 
}
\item{ibd.var }{
a list containing an element for each pedigree. The elements include
ped.id, person1.id, person2.id, and the variance-covariance matrix 
of ibd statistics between pairs of subject pairs.
}
\item{status}{
column name of data that has affection status (1=unaffected;
2=affected; NA=missing)
}
\item{ped.id }{
column name of data that has pedigree id
}
\item{person.id }{
column name of data that has a person's id, which only have to be unique
within pedigrees
}
\item{data}{
a data frame containing all variables in the model formula, in addition to 
status, ped.id, and person.id
}
\item{c.scale  }{
choice for scaling factor on covariates, "minimax" (same as LODPAL) 
and "nodom" (no dominance-variance)
}
\item{subset}{
an optional vector specifying a subset of observations in data
}
\item{weights}{
an optional vector of weights to be used in the regression method
}
\item{na.action }{
a function that defines how missing values (NA) are handled.  
The default is set by the 'na.action' setting of 'options', 
and is 'na.fail' by default.
}
\item{min.pairs }{
minimum number of relative pairs needed for linkage tests to be
calculated on a specific status group (e.g. AA pairs)
}
\item{epsilon}{
minimum value for singular values in generalized inverse calculations
}
\item{\dots }{
further arguments passed to or from other methods
}
}
\value{
A list with the following components:

\item{Call}{
function call to ibdreg 
}
\item{relpair.tbl}{
table with number of relative pairs in each pedigree
}
\item{status.tbl}{
table of counts for AA, UU, and AU relative pairs 
}
\item{AA.linkage}{
an object of class linkage.tests; includes score statistics for linkage 
for AA relative pairs that includes tests (i) - (iii), (v) (see details)
}
\item{UU.linkage}{
an object of class linkage.tests; includes score statistics for linkage 
for UU relative pairs that includes tests (i) - (iii), (v) (see details)
}
\item{AU.linkage}{
an object of class linkage.tests; includes score statistics for linkage 
for AU relative pairs that includes tests (i) - (iii), (v) (see details)
}
\item{ALL.linkage}{
object of class linkage.all; includes score statistics for linkage
with ALL pairs (test iv).
}
}
\details{
The tests for genetic linkage use quasi-likelihood score statistics, 
formulated in terms of weighted least squares regression.  
The covariates in the regression framework are scaled according to 
the degree of relationship between relative pairs within the pedigree 
(c.scale). The method yields the following tests for linkage and/or 
covariate effect for relative pairs:  (i) Linkage only;  
(ii) Linkage with covariate effects; (iii) Covariate effect 
on IBD sharing (e.g, heterogeneity), assuming either a 
model-based variance-covariance matrix or a robust variance-covariance
matrix. Tests (i) - (iii) are evaluated as unconstrained (two-sided test) and 
constrained for excess allele sharing in one direction.  If formula is
\~1, there are no covariates, so the tests will be for linkage only.

Another test (iv) is for linkage using all relative pairs. This test 
imposes constraints on allele sharing for one-sides tests favoring
linkage, achieved by constraining IBD allele sharing for pairs as
AA > UU > null, and AU < null, where null is the expected IBD sharing
without linkage.

Tests (i) - (iv) are performed using quasi-likelihood scores, 
and their distributions are asymptotically chi-square or mixture chi-square.

Another set of tests (v) are z-scores measuring the departure from
expected allele sharing for each pedigree, at each position.  These
are provided for diagnostic purposes in the linkage.tests results for
AA, UU, and AU pairs.

Correlations of relative pairs from the same pedigree are accounted for by a 
covariance matrix for ibd statistics within pedigrees.  This
variance-covariance matrix is calculated using results from merlin
(with options "--ibd --matrices").
}
\section{References}{
Schaid DJ, Sinnwell JP, Thibodeau SN. Testing Genetic
Linkage with Affected Relative Pairs and Covariates by
Quasi-Likelihood Score Statistics. Submitted.
}
\seealso{
\code{\link{create.ibd.dat}},
\code{\link{sim.ibd.var}},
\code{\link{exact.ibd.var}},
\code{\link{linkage.tests}},
\code{\link{linkage.all}}, 
\code{\link{print.ibdreg}}, 
\code{\link{plot.ibdreg}}
}
\examples{
## see manual for examples, the data requirements are here:

# make a data.frame with ped.id, person.id, status, and covariates
# make ibd.dat object using create.ibd.dat
# make ibd.var object using create.ibd.var

# call ibdreg
}
% docclass is function
% Converted by Sd2Rd version 43267.
