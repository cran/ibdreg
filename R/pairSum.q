#$Author: sinnwell $
#$Date: 2006/09/12 21:25:23 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/pairSum.q,v 1.1 2006/09/12 21:25:23 sinnwell Exp $
#$Locker:  $
#$Log: pairSum.q,v $
#Revision 1.1  2006/09/12 21:25:23  sinnwell
#Initial revision
#
#Revision 1.1  2006/04/20 18:32:11  sinnwell
#Initial revision
#

## Jason Sinnwell   4/2006
## Mayo Clinic, Division of Biostatistics
## project: ibdreg software for SPLUS/R

pairSum <- function(cov1,cov2) {
  return(cov1 + cov2)
}
