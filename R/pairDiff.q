#$Author: sinnwell $
#$Date: 2006/09/12 21:25:33 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/pairDiff.q,v 1.1 2006/09/12 21:25:33 sinnwell Exp $
#$Locker:  $
#$Log: pairDiff.q,v $
#Revision 1.1  2006/09/12 21:25:33  sinnwell
#Initial revision
#
#Revision 1.1  2006/04/20 18:32:03  sinnwell
#Initial revision
#

## Jason Sinnwell   4/2006
## Mayo Clinic, Division of Biostatistics
## project: ibdreg software for SPLUS/R

pairDiff <- function(cov1,cov2) {

  return( abs(cov1 - cov2) )

}

