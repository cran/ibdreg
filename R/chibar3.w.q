#$Author: sinnwell $
#$Date: 2006/03/08 16:43:20 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/chibar3.w.q,v 1.1 2006/03/08 16:43:20 sinnwell Exp $
#$Locker:  $
#$Log: chibar3.w.q,v $
#Revision 1.1  2006/03/08 16:43:20  sinnwell
#Initial revision
#

# Dan Schaid
# November, 2005
# Mayo Clinic, Division of Biostatistics

chibar3.w <- function(sigma){

  # compute mixture proportions, w, for mixture chi-square
  # distribution when k=0..3, where k = df of chi-squares
  # sigma = covariance matrix (3 x 3)
  # elements of returned vector w are mixing proportions, with
  # element  w[k+1] the mixture proportion for k, k=0..3

  se <- sqrt(diag(sigma))
  r <- sigma/(se %o% se)

  r.partial <- function(r, i, j, k){
    r.partial <- (r[i,j] - r[i,k]*r[j,k])/(sqrt(1-r[i,k]^2) * sqrt(1-r[j,k]^2))
    return(r.partial) 
  }

  w <- numeric(4)

  w[4] <- (2*pi - acos(r[1,2]) - acos(r[1,3]) - acos(r[2,3]))/(4*pi)
  w[3] <- (3*pi - acos(r.partial(r,1,2,3)) - acos(r.partial(r,1,3,2)) -
           acos(r.partial(r,2,3,1)))/(4*pi)
  w[2] <- .5 - w[4]
  w[1] <- .5 - w[3]

  return(w)

}
