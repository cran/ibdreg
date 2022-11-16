% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:39:33 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.ibd.var.Rd,v 1.2 2010/02/08 sinnwell Exp $ -->

\name{sim.ibd.var}
\alias{sim.ibd.var}
\title{
  Create an ibd.var object via simulations using gene-dropping
}
\description{
Create an ibd.var object, containing a mean vector and covariance
matrix for the ibd sharing values between each pair of relative
pairs in each pedigree.  Use a gene-dropping approach with 2 unique
alleles for each founder, and simulate n.sim times.  Then calculate the
moments based on the vector of allele sharing between all relative 
pairs within the pedigree.
}
\usage{
sim.ibd.var(pedfile, 
            male.code=1, 
            female.code=2, 
            x.linked=FALSE, 
            n.sim=1000, 
            print.resources=FALSE)
}
\arguments{

\item{pedfile}{
file containing ids for pedigree, person, father, and mother, and
gender.  Also known as a pre-MAKEPED format file.
}
\item{male.code }{
The code that denotes a sex classification of "male".
}
\item{female.code }{
The code that denotes a sex classification of "female".
}
\item{x.linked }{
A flag that denotes whether the inheritance mode should be autosomal
(x.linked=FALSE, the default) or X-linked (x.linked=TRUE).
}
\item{n.sim }{
Number of simulations
}
\item{print.resources }{
logical, set to TRUE to print CPU and memory usage for
calculations on each pedigree.
}
}
\value{
The return value, call it ret, is an ibd.var object that contains 
the following elements for each pedigree:
\item{ped.id}{ pedigree id}
\item{person1.id}{ vector of ids for first person in the relative pair}
\item{person2.id}{ vector of ids for second person in the relative pair}
\item{sm}{ mean vector of ibd sharing between relative pairs:  person1.id, person2.id}
\item{sv}{ variance-covariance matrix for ibd sharing between pairs of relative pairs}
}
\details{
The bulk of the work for simulating ibd sharing vectors is in
sim.mark.prop(). These methods will handle practically
any correctly formed pedigree structure.  This includes pedigrees 
containing loops, multiple pairs of ancestoral founders, and both 
loops and multiple pairs of ancestoral founders.  Thus, maximum 
flexibility may be obtained.  Both autosomal and X-linked methods 
of inheritance are supported through the use of the x.linked parameter. 
}
\section{References}{
Schaid DJ, Sinnwell JP, Thibodeau SN. Testing Genetic Linkage with
Relative Pairs and Covariates by Quasi-Likelihood Score Statistics.  Submitted.
}
\seealso{
\code{\link{exact.ibd.var}},
\code{\link{sim.ibd.setup}},
\code{\link{sim.mark.prop}}
}
\examples{
## See manual(s) for full example usage
## Below is a basic example for one pedigree
## Not Run:
# Since input parameter is a "pre" file, create a "pre" file for one pedigree

# make a data.frame to write to a file
ped.id <- rep(7,7)
person <- c(1,2,3,4,5,6,7)
father <- c(0,0,0,1,0,3,5)
mother <- c(0,0,0,2,0,4,4)
sex    <- c(1,2,1,2,1,1,1)
chrom1 <- c(1,3,5,0,7,0,0)
chrom2 <- c(2,4,6,0,8,0,0)
ped7.df <- data.frame(ped.id, person, father, mother, 
                     sex, chrom1, chrom2)

# write the file
\dontrun{write.table(ped7.df, file="ped7.pre", row.names=FALSE, col.names=FALSE)}

\dontrun{ped7.ibdVar <- sim.ibd.var("ped7.pre", n.sim=1000)}

# results are long and difficult to interpret
# here is an example of how to view
\dontrun{print(ped7.ibdVar[[1]], digits=2)}

# demonstrate for x.linked=TRUE
\dontrun{ped7.ibdVarX <- sim.ibd.var("ped7.pre", n.sim=1000, x.linked=TRUE)}
}
% docclass is function
% Converted by Sd2Rd version 43267.
