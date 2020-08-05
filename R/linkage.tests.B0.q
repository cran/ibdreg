#$Author: sinnwell $
#$Date: 2006/03/08 16:37:23 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.tests.B0.q,v 1.1 2006/03/08 16:37:23 sinnwell Exp $
#$Locker:  $
#$Log: linkage.tests.B0.q,v $
#Revision 1.1  2006/03/08 16:37:23  sinnwell
#Initial revision
#
######################################
#   Jason Sinnwell,  Daniel Schaid
#   Div of Biostatistics
#   Mayo Clinic, HSR  2005
######################################

linkage.tests.B0 <- function(y.mat, xvec, ibdvar.lst) {
  
  # pre-calculate the mean ibd sharing, the intercept for the linear model

  cvy.mat <- matrix(0,nrow=length(ibdvar.lst), ncol=ncol(y.mat))
  cvc.vec <- rep(0, ncol(y.mat))
  ped.start <- 0

  for(ped in 1:length(ibdvar.lst)) {

     # row number corresponding to persons in pedigree ped
    this.ped <- (ped.start+1):(ped.start+ibdvar.lst[[ped]]$n)
    if(!is.na(ibdvar.lst[[ped]]$rank) & length(this.ped)>0) {
          #cat('ped', ped, '\n')

      # save portions of Bo hat calculation, contributed by each ped
      cvy.mat[ped,] <- t(xvec[this.ped])%*%ibdvar.lst[[ped]]$sv.ginv%*%y.mat[this.ped,,drop=FALSE]
      cvc.vec[ped] <- t(xvec[this.ped])%*%ibdvar.lst[[ped]]$sv.ginv%*%xvec[this.ped]
    }
    # re-set ped.start to be last person from ped
    ped.start <- ped.start+ibdvar.lst[[ped]]$n

  }
 
  apply(cvy.mat, 2, sum)/sum(cvc.vec, na.rm=TRUE)

}
