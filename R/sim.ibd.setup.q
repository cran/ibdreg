#$Author: sinnwell $
#$Date: 2006/10/10 21:20:38 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.ibd.setup.q,v 1.1 2006/10/10 21:20:38 sinnwell Exp $
#$Locker:  $
#$Log: sim.ibd.setup.q,v $
#Revision 1.1  2006/10/10 21:20:38  sinnwell
#Initial revision
#
## authors: Dan Schaid and Dan Folie
## purpose: within sim.ibd.var, setup up a ped object for sim.mark.prop
## package: ibdreg
## date:    Oct 2006

sim.ibd.setup <- function(ped, miss.val = c(NA, 0), x.linked=FALSE)
{

	n <- length(ped$person)

        a1 <- rep(1,n)
        a1.orig <- rep(1,n)
        a2 <- rep(1,n)
        a2.orig <- rep(2,n)

        ped <- cbind(ped, a1=a1, a1.orig=a1.orig, a2=a2, a2.orig=a2.orig)

        ped <- setFounderAlleles(ped, x.linked=x.linked)

        return(ped)
}
