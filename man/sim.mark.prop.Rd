% $Author: sinnwell $ -->
% $Date: 2010/02/08 18:39:33 $ -->
% $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.mark.prop.Rd,v 1.2 2010/02/08 18:39:33 sinnwell Exp $ -->
% $Locker:  $ -->
% 
%$Log: sim.mark.prop.Rd,v $
%Revision 1.2  2010/02/08 18:39:33  sinnwell
%rm empty sections
%
%Revision 1.1  2008/10/14 18:09:51  sinnwell
%Initial revision
%
%Revision 1.6  2006/10/30 21:40:36  sinnwell
%remove pre tags
%
%Revision 1.5  2006/10/30 21:14:29  sinnwell
%reduce example
%
%Revision 1.4  2005/10/07 18:23:16  folie
%Updated SGML file
%
%Revision 1.3  2005/10/04 18:42:48  folie
%Edited comments
%
%Revision 1.2  2005/10/03 21:25:57  folie
%Updated SGML help file
%
%Revision 1.1  2005/09/28 18:37:57  folie
%Initial revision
%
\name{sim.mark.prop}
\alias{sim.mark.prop}
\title{
  Simulate Marker Propagation
}
\description{
Simulates inheritance-by-descent in pedigree data through both autosomal
and X-linked modes of inheritance.
}
\usage{
sim.mark.prop(ped, 
              iseed=NULL, 
              miss.val=c(NA, 0), 
              male.code=1, 
              female.code=2, 
              x.linked=FALSE, 
              proband=0, 
              n.iter=1)
}
\arguments{
\item{ped}{
The ped argument is an object of class "list" or class "data.frame" and
it must have the following six fields:

ped\$person:  A vector or factor.  The person ID of each pedigree member.

ped\$father:  A vector or factor.  The father ID of each pedigree member.

ped\$mother:  A vector or factor.  The mother ID of each pedigree member.

ped\$sex:     A vector or factor.  The sex of each pedigree member.

ped\$chrom1:  A vector or factor.  The first marker of the marker data 
             pair for each pedigree member (only founder data is required;
             non-founder data will be ignored in the simulations).  If 
             x.linked=TRUE, then male pedigree members must be homozygous
             at the specified marker site.

ped\$chrom2:  A vector or factor.  The second marker of the marker data
             pair for each pedigree member (only founder data is required;
             non-founder data will be ignored in the simulations).  If 
             x.linked=TRUE, then male pedigree members must be homozygous
             at the specified marker site.
}
\item{iseed}{
An integer or a saved copy of .Random.seed.  This allows simulations to be
reproduced by using the same initial seed.  If no value is given, the
current copy of .Random.seed will be used.
}
\item{miss.val }{
Codes that denote missing values in the pedigree data given by ped.
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
\item{proband}{
An integer denoting which person should be the first person encountered
by the traversal algorithm.  This parameter is used for testing and need
not be set by the user.  It's value will not have any effect on the 
statistical properties of the return value.
}
\item{n.iter }{
An integer stipulating the number of simulated inheritance-by-descent
calculations that should be computed (and returned) by the 
sim.mark.prop() program.
}
}
\value{
The return value, call it ret, is an object of class "list" or class
"data.frame" (the same class as the class of the "ped" parameter that
the sim.mark.prop() function was called with, and it has the follwing
fields:

ret\$person:   A vector or factor.  The person ID of each pedigree member.

ret\$father:   A vector or factor.  The father ID of each pedigree member.

ret\$mother:   A vector or factor.  The mother ID of each pedigree member.

ret\$sex:      A vector or factor.  The sex of each pedigree member.

ret\$chrom1:   A vector or factor.  The first marker of the marker data 
              pair for each pedigree member input into the 
              sim.mark.prop() program.

ret\$chrom2:   A vector or factor.  The second marker of the marker data 
              pair for each pedigree member input into the 
              sim.mark.prop() program.

ret\$mark1:    A vector if n.iter=1, otherwise an object of class 
              "model.matrix".  Letting n.subjects denote the number of
              pedigree members, "mark1" is an n.subject-by-n.iter matrix
              with the ith column listing the first simulated marker of
              each simulated marker pair computed by the ith simulation
              for the corresponding pedigree member.  If X-linked 
              inheritance was specified (x.linked=TRUE), each male will
              be homozygous at the simulated marker location.

ret\$mark2:    A vector if n.iter=1, otherwise an object of class 
              "model.matrix".  Letting n.subjects denote the number of
              pedigree members, "mark2" is an n.subject-by-n.iter matrix
              with the ith column listing the second simulated marker of
              each simulated marker pair computed by the ith simulation
              for the corresponding pedigree member.  If X-linked 
              inheritance was specified (x.linked=TRUE), each male will
              be homozygous at the simulated marker location

ret\$x.linked: A flag denoting whether the simulated "mark1" and "mark2"
              data was simulated to conform to autosomal inheritance or
              X-linked inheritance.
}
\section{Side Effects}{
The sim.mark.prop() program has no side effects.
}
\details{
The sim.mark.prop() program will handle practically any correctly formed
pedigree structure.  This includes pedigrees containing loops, multiple 
pairs of ancestoral founders, and both loops and multiple pairs of 
ancestoral founders.  Thus, maximum flexibility may be obtained.

Both autosomal and X-linked methods of inheritance are supported through
the use of the x.linked parameter.  When X-linked inheritance is 
specified, the input marker data ("ped\$chrom1" and "ped\$chrom2") must be
homozygous for male pedigree members, and the simulated output marker 
data ("ret\$mark1" and "ret\$mark2") will be homozygous for male pedigree
members.

When more than one set of simulated marker data is desired, it is
advisable to set the "n.iter" parameter to the number of simulated
marker sets that are desired, rather than calling sim.mark.prop() in a
loop n.iter times, as the underlying compiled C program that does the 
simulations forms the inputted pedigree one time, and then perfoms 
n.iter simulations upon that single pedigree.  Calling sim.mark.prop()
one time for each simulation results in proportionately large increases
in redundant work (both in S-PLUS and C) and call overhead.  The only 
thing to be aware of is that the amount of memory used by a call to 
sim.mark.prop() grows linearly with (n.subjects*n.iter).  Therefore, if
the product of the number of pedigree subjects with the number of 
simulations requested is unusually large, it may be advisable to split
the simulations into more than one call to sim.mark.prop() in order to 
avoid exhausting the amount of primary memory afforded to S-PLUS (both
to avoid the unintended termination of the host S-PLUS process and to 
avoid the significant increase in computational time that is seen when
disk paging is required).

The non-founder members of the pedigree will have the first marker of 
each simulated marker pair (ret\$mark1) inherited from the father, and
the second marker of each simulated marker pair (ret\$mark2) inherited 
from the mother.  Thus, it is possible to follow maternal and paternal
allele transmission through dependent pedigree members.
}
\examples{
###################################################
### The following examples use a list as input: ###
###################################################

ped7 <- list(NULL)
ped7$person <- c(1,2,3,4,5,6,7)
ped7$father <- c(0,0,0,1,0,3,5)
ped7$mother <- c(0,0,0,2,0,4,4)
ped7$sex    <- c(1,2,1,2,1,1,1)
ped7$chrom1 <- c(1,3,5,0,7,0,0)
ped7$chrom2 <- c(2,4,6,0,8,0,0)
 
#ped7:

#       1-------2       
#     1/2   |   3/4     
#           |           
#           |           
#   3-------4-------5   
# 5/6   |  0/0  |   7/8 
#       |       |       
#       |       |       
#       6       7       
#     0/0       0/0     

prop7 <- sim.mark.prop(ped7, n.iter=5)

prop7


## an x.linked example

ped7.x <- list(NULL)
ped7.x$person <- c(1,2,3,4,5,6,7)
ped7.x$father <- c(0,0,0,1,0,3,5)
ped7.x$mother <- c(0,0,0,2,0,4,4)
ped7.x$sex    <- c(1,2,1,2,1,1,1)
ped7.x$chrom1 <- c(1,2,4,0,5,0,0)
ped7.x$chrom2 <- c(1,3,4,0,5,0,0)

# ped7.x:
#
#       1-------2       
#     1/1   |   2/3     
#           |           
#           |           
#   3-------4-------5   
# 4/4   |  0/0  |   5/5 
#       |       |       
#       |       |       
#       6       7       
#     0/0       0/0 

prop7X <- sim.mark.prop(ped7.x, x.linked=TRUE, n.iter=5)

prop7X
}
% docclass is function
% Converted by Sd2Rd version 43267.
