#$Author: sinnwell $
#$Date: 2006/11/07 21:34:29 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/create.pairs.frame.cvec.q,v 1.3 2006/11/07 21:34:29 sinnwell Exp $
#$Locker:  $
#$Log: create.pairs.frame.cvec.q,v $
#Revision 1.3  2006/11/07 21:34:29  sinnwell
#no.dom to nodom
#
#Revision 1.2  2006/03/15 20:33:44  sinnwell
#a1, a2 to prior1, prior2
#
#Revision 1.1  2006/03/08 16:46:25  sinnwell
#Initial revision
#

## Jason Sinnwell and Dan Schaid
## 11/2005
## Mayo Clinic Division of Biostatistics


create.pairs.frame.cvec <- function(c.scale, prior2, prior1) {
  # internal use
  # two choices for c.scale, minimax (lodpal's) or nodom (no-dominance)

  if( c.scale == "minimax" ) {
    c.vec <-  7.268*prior2 - 5.634*prior2*prior1 -
              7.268*prior2^2 + prior1 - prior1^2
  } else if (c.scale=="nodom") {
    mu <- 2*prior2 + prior1
    c.vec <- 2*prior2 * (2-mu) + prior1*(1-mu)
  }
  c.vec
}
