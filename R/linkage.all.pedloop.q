#$Author: sinnwell $
#$Date: 2006/05/17 19:55:17 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.all.pedloop.q,v 1.2 2006/05/17 19:55:17 sinnwell Exp $
#$Locker:  $
#$Log: linkage.all.pedloop.q,v $
#Revision 1.2  2006/05/17 19:55:17  sinnwell
#take out unexpected stats
#
#Revision 1.1  2006/03/08 16:37:23  sinnwell
#Initial revision
#

######################################
#   Jason Sinnwell,  Daniel Schaid
#   Div of Biostatistics
#   Mayo Clinic, HSR  2005
######################################

linkage.all.pedloop <- function(y.vec, x.adj, ibdvar.lst) {
  
  # LOOP OVER PEDIGREES, computing linkage test stats specific to 1 position
  # save contributions of u.score from pedigrees to linkage.all
  # for further checks on unexpected allele-sharing

  p <- ncol(x.adj)

  ##---- INITIALIZE DATA OBJECTS------
  # add to them throughout the ped loop
  # keep contributions of u.score in u.score.ped
  # "                  "  u.var in u.var.ped (diagonal only--variances)
  u.var <- matrix(0,ncol=p,nrow=p)
  u.score <- matrix(0,ncol=p,nrow=1)
  u.score.ped <- matrix(0,ncol=p, nrow=length(ibdvar.lst))
  ##  u.var.ped <- matrix(0,ncol=p, nrow=length(ibdvar.lst))
  
  ##-------ITERATE OVER PEDIGREES--------

  ped.start <- 0
  for(i in 1:length(ibdvar.lst)) {

    # create index for current pedigree
    this.ped <- (ped.start+1):(ped.start+ibdvar.lst[[i]]$n)

    # if no rank of the ibdvar mat, can't do anything
    if(!is.na(ibdvar.lst[[i]]$rank)) {

      y.sub <- y.vec[this.ped]
      x.sub <- x.adj[this.ped,,drop=FALSE] 
      Vo.ginv <- ibdvar.lst[[i]]$sv.ginv

      u.score.ped[i,] <- t(y.sub)%*%Vo.ginv%*%x.sub
#      u.var.ped[i,] <- diag(t(x.sub)%*%Vo.ginv%*%x.sub)

      # sum score and variance for first test, sum over pedigrees
      u.score <- u.score + u.score.ped[i,]    # t(y.sub)%*%Vo.ginv%*%x.sub
      u.var <- u.var + t(x.sub)%*%Vo.ginv%*%x.sub

    }  # if !is.na(rank)

    ped.start <- ped.start + ibdvar.lst[[i]]$n

  }

  # RETURN OBJECT
  obj <- list(u.score = u.score,
              u.var = u.var)

  return(obj)
}
