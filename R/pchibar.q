#$Author: sinnwell $
#$Date: 2006/03/08 16:45:25 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/pchibar.q,v 1.1 2006/03/08 16:45:25 sinnwell Exp $
#$Locker:  $
#$Log: pchibar.q,v $
#Revision 1.1  2006/03/08 16:45:25  sinnwell
#Initial revision
#
# Dan Schaid
# November, 2005
# Mayo Clinic, Division of Biostatistics

pchibar <- function(x, df, wt){

  # compute P(X <= x), where P is the chi-bar distribution
  # df = vector of df
  # wt = vector of mixing proportions

  if(x<=0){
      return(0)
  }

  zed <- df==0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pchisq(x, df[!zed])*wt[!zed])

  return(cdf)
}
