\name{create.ibd.dat}
\alias{create.ibd.dat}
\title{
  Create an ibd.dat object
}
\description{
Combine a posterior IBD probabilities file with a prior
IBD probabilities file
}
\usage{
create.ibd.dat(postfile, priorfile, software="merlin", x.linked=FALSE,
               cov.data=NULL, rm.noninform=TRUE)
}
\arguments{
\item{postfile}{
full path and name of the file with posterior IBD probabilities
}
\item{priorfile}{
full path and name of the file with prior IBD probabilities
}
\item{software}{
character string of which software was used to create IBD probability files
}
\item{x.linked }{
logical, is the chromosome X-linked.  If TRUE, cov.data is required
}
\item{cov.data }{
name of data.frame containing covariates, specifically  ped.id,
person.id, and sex.  Required when x.linked=TRUE.
}
\item{rm.noninform }{
logical, if TRUE, remove relative pairs that are not informative for linkage
}
}
\value{
An object with class ibd.dat, which contains the following
elements:

\item{ped.id}{
pedigree identifier code
}
\item{person1.id}{
identifier to person 1 of the relative pair
}
\item{person2.id}{
identifier to person 2 of the relative pair
}
\item{post0}{
data.frame with probability of sharing zero (0) alleles ibd
between relative pairs (rows) at each position (columns)
}
\item{post1}{
data.frame with probability of sharing one (1) alleles ibd
between relative pairs (rows) at each position (columns)
}
\item{post2}{
data.frame with probability of sharing two (2) alleles ibd
between relative pairs (rows) at each position (columns)
}
\item{prior0}{
vector with probability of sharing zero (0) alleles ibd
between relative pairs, given no genotype data
}
\item{prior1}{
vector with probability of sharing one (1) alleles ibd
between relative pairs, given no genotype data
}
\item{prior2}{
vector with probability of sharing two (2) alleles ibd
between relative pairs, given no genotype data
}
}
\section{Side Effects}{

}
\details{
Perl scripts are provided within ibdreg/perl/ for creating IBD
probability files, and are explained in the user manual.  
Prior probability files are made using a homozygous marker 
on the pedigree structure, or a "dummy" marker.  
When the chromosome is X-linked, the Merlin software treats males as
homozygous for their X chromosome, and thus the probability of sharing 1
and 2 alleles IBD with any other relative are switched within 
create.ibd.dat.
}
\section{References}{

}
\seealso{
\code{\link{ibdreg}},
\code{\link{ibd.df.merlin}},
\code{\link{mergeIBD}}
}
\examples{
## do not run example in testing
## uncomment to run for demo

## ibd file for 1 chromosome
# ibdfile.ch20 <- "post.ibd"  

## ibd file for 1 locus for prior probs.
# prior.ibdfile <- "prior.ibd"

# ibd.dat.obj <- create.ibd.dat(postfile=post.ibd, 
#                               priorfile=prior.ibd,
#                               rm.noninform=FALSE)

# names(ibd.dat.obj)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
