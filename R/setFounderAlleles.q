#$Author: sinnwell $
#$Date: 2006/10/13 16:44:43 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/setFounderAlleles.q,v 1.2 2006/10/13 16:44:43 sinnwell Exp $
#$Locker:  $
#$Log: setFounderAlleles.q,v $
#Revision 1.2  2006/10/13 16:44:43  sinnwell
#ped$male did not exist, but ped$sex did.
#
#Revision 1.1  2006/10/10 21:21:47  sinnwell
#Initial revision
#
## authors: Dan Schaid and Dan Folie
## purpose: set founder alleles for a pedigree
## package: ibdreg
## date:    Oct 2006

setFounderAlleles <- function(ped, x.linked=FALSE)
{

# changed chrom1 and chrom2 default missing values to 0 (rMarkFounders has NA)

     #identify founders
	founder <- ped$father == 0 & ped$mother == 0
        male <- ped$sex==1

	# n = number of alleles

        n.male <- sum(founder & male)
        n.female <- sum(founder & (!male) )


	tmp1 <- tmp2 <- NULL
	nSubj <- length(founder)
	nMark <- 1

        if(x.linked){

            n.allele <- n.male + 2*n.female

            a1.male <- 1:n.male
            a2.male <- 1:n.male

            a1.female <- seq((n.male+1), n.allele, by=2)
            a2.female <- seq((n.male+2), n.allele, by=2)

         } else {

            n.allele <- 2*(n.male + n.female)

            a1.male <- seq(1,2*n.male, by=2)
            a2.male <- seq(2,2*n.male, by=2)

            a1.female <- seq((2*n.male+1), n.allele, by=2)
            a2.female <- seq((2*n.male+2), n.allele, by=2)
         }

	chrom1 <- chrom2 <- matrix(rep(0, nSubj * nMark), nrow = nSubj)

	chrom1[founder & male,  ] <- a1.male
        chrom1[founder & (!male), ] <- a1.female
	sr.class(chrom1) <- "model.matrix"

	chrom2[founder & male,  ] <- a2.male
        chrom2[founder & (!male), ] <- a2.female

	sr.class(chrom2) <- "model.matrix"
	ped <- data.frame(ped, chrom1 = chrom1, chrom2 = chrom2)
	return(ped)
}
