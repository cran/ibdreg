#$Author: sinnwell $
#$Date: 2006/05/17 20:22:33 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.tests.pedloop.q,v 1.3 2006/05/17 20:22:33 sinnwell Exp $
#$Locker:  $
#$Log: linkage.tests.pedloop.q,v $
#Revision 1.3  2006/05/17 20:22:33  sinnwell
#include z.ped for ped sharing departures
#
#Revision 1.2  2006/04/26 14:13:34  sinnwell
#add epsilon
#
#Revision 1.1  2006/03/08 16:37:23  sinnwell
#Initial revision
#

######################################
#   Jason Sinnwell,  Daniel Schaid   
#   Div of Biostatistics
#   Mayo Clinic, HSR  2005
######################################

linkage.tests.pedloop <- function(y.vec, x.adj, ibdvar.lst, B0, concordant=TRUE, epsilon) {

  # LOOP OVER PEDIGREES, computing linkage test stats specific to 1 position
  # concordant pairs: AA or UU, both restricted B0>0
  # discordant pairs: AU restricted B0<0
  
    ## variables needed:
    ##  Bo.hat[pos], y.mat[,pos], x.adj, ibdvar.lst, pos
    ##
    ## cons == constrained test, constrain B0 or u.score[1] > 0

  ncov <- ncol(x.adj)

 ##----- INITIALIZE ALL DATA OBJECTS------
 ##  add to them throughout the ped loop
  u.var <- matrix(0,ncol=ncov,nrow=ncov)
  u.score <- matrix(0,ncol=ncov,nrow=1)
  z.ped <- numeric(length(ibdvar.lst))
                         
  if(ncov>1) {
    u.noint <- matrix(0, nrow=1, ncol=ncov-1)
    u.noint.cons <- matrix(0, nrow=1, ncol=ncov-1)
    u.int <- matrix(0, nrow=1,ncol=1)
    u.int.cons <- matrix(0, nrow=1,ncol=1)
  
    sum.xvx <- matrix(0,ncol=ncov-1,nrow=ncov-1)
    sum.xvc <- matrix(0,nrow=ncov-1,ncol=1)
    sum.cvc <- 0
    
    sum.u2u2 <- matrix(0,nrow=ncov-1,ncol=ncov-1)
    sum.u2u1 <- matrix(0,nrow=ncov-1,ncol=1)
    sum.u1u1 <- matrix(0,ncol=1,nrow=1)
    sum.u2u2.cons <- matrix(0,nrow=ncov-1,ncol=ncov-1)
    sum.u2u1.cons <- matrix(0,nrow=ncov-1,ncol=1)
    sum.u1u1.cons <- matrix(0,ncol=1,nrow=1)

    if(concordant) B0.cons <- ifelse(B0<0, 0, B0)
    else B0.cons <- ifelse(B0>0, 0, B0)
    
  } # if ncov>1
  
  #-------ITERATE OVER PEDIGREES--------

  this.ped.start <- 0
  for(i in 1:length(ibdvar.lst)) {

    # create index for current pedigree
    this.ped <- (this.ped.start+1):(this.ped.start+ibdvar.lst[[i]]$n)

    # if no rank of the ibdvar mat, can't do anything
    if(!is.na(ibdvar.lst[[i]]$rank)) {

      y.sub <- y.vec[this.ped]
      x.sub <- x.adj[this.ped,,drop=FALSE] 
      Vo.ginv <- ibdvar.lst[[i]]$sv.ginv

      # sum score and variance for first test, sum over pedigrees
      u.score.ped <- t(y.sub) %*% Vo.ginv %*% x.sub
      u.var.ped <- t(x.sub) %*% Vo.ginv %*% x.sub
      u.score <- u.score + u.score.ped
      u.var <- u.var + u.var.ped

      z.ped[i] <- u.score.ped[1]/sqrt(u.var.ped[1,1])
      
      if(ncov>1) {
       # divide the score vector:
       #   noint => without B0, int => intercept only
        
        # keep the ith portion for below, and cumulative sum
        u.noint.i <- t(y.sub-B0*x.sub[,1])%*%Vo.ginv%*%x.sub[,-1, drop=FALSE]
        u.noint.i.cons <-  t(y.sub-B0.cons*x.sub[,1])%*%Vo.ginv%*%x.sub[,-1, drop=FALSE]

        u.noint <- u.noint + u.noint.i
        u.noint.cons <- u.noint.cons + u.noint.i.cons 

        u.int.i <- t(y.sub-B0*x.sub[,1])%*%Vo.ginv%*%x.sub[,1] 
        u.int.i.cons <- t(y.sub-B0.cons*x.sub[,1])%*%Vo.ginv%*%x.sub[,1] 
        u.int <- u.int + u.int.i
        u.int.cons <- u.int.cons + u.int.i.cons

       # model-based variance matrix has 3 parts
        sum.xvx <- sum.xvx + t(x.sub[,-1, drop=FALSE])%*%Vo.ginv%*%x.sub[,-1, drop=FALSE]
        sum.xvc <- sum.xvc + t(x.sub[,-1, drop=FALSE])%*%Vo.ginv%*%x.sub[,1]
        sum.cvc <- sum.cvc + t(x.sub[,1])%*%Vo.ginv%*%x.sub[,1]

       # robust variance matrix has 3 parts, partitioned into u2u2, u2u1=u1u2, u1u1
        sum.u2u2 <- sum.u2u2 + t(u.noint.i)%*%u.noint.i
        sum.u2u1 <- sum.u2u1 + t(u.noint.i)%*%u.int.i
        sum.u1u1 <- sum.u1u1 + u.int.i*u.int.i

      # constrained
        sum.u2u2.cons <- sum.u2u2.cons + t(u.noint.i.cons)%*%u.noint.i.cons
        sum.u2u1.cons <- sum.u2u1.cons + t(u.noint.i.cons)%*%u.int.i.cons
        sum.u1u1.cons <- sum.u1u1.cons + u.int.i.cons*u.int.i.cons

      } # if ncov>1

    }  #if !is.na(rank)

    this.ped.start <- this.ped.start + ibdvar.lst[[i]]$n

  }

  ## LINKAGE ONLY, NO COVARIATE
     # not constrained; d.f. assumed to be 1
  if(u.var[1,1]==0) stop("Variance of zero(0) in denominator")
  test.linkage <- u.score[1]*u.score[1]/u.var[1,1]  #t(u)%*%inv(v)%*%u

    # constrained, a one-sided test, mixture distribution
    # d.f. assumed to be 50/50 mix of (0,1)
  if(concordant) test.linkage.cons <- ifelse(u.score[1]>0, test.linkage, 0)
  else test.linkage.cons <- ifelse(u.score[1]<0, test.linkage, 0)
  
  if(ncov>1) {  
    ############ LINKAGE WITH COVARIATE ################
    # not constrained
    uvar.Ginv <- Ginv(u.var, eps=epsilon)
    dof.linkcov <- dof2.linkcov.cons <- uvar.Ginv$rank
    dof1.linkcov.cons <- dof2.linkcov.cons-1
    test.linkcov <- u.score%*%uvar.Ginv$Ginv%*%t(u.score)
    
    # constrained statistics
    if((concordant & u.score[1] < 0) || (!concordant & u.score[1]>=0)) {
      u.var.00 <- u.var[1,1] # checked for eq 0 above
      u.var.10 <- u.var[-1,1, drop=FALSE]
      u.score.cons <- u.score[-1] - u.var.10 * (u.score[1]/u.var.00)
      u.var.cons.Ginv <- Ginv(u.var[-1,-1, drop=FALSE], eps=epsilon)

      test.linkcov.cons <- t(u.score.cons)%*%u.var.cons.Ginv$Ginv%*%u.score.cons
      dof1.linkcov.cons <- u.var.cons.Ginv$rank  # override above dof1, in this case
      
    } else {
      test.linkcov.cons <- test.linkcov
    }

    ###### COVARIATE EFFECT TEST WITH MODEL-BASED VARIANCE #####
    # not constrained
    u.var.mod <- sum.xvx - (sum.xvc)%*%(1/sum.cvc)%*%t(sum.xvc )
    u.var.mod.ginv <- Ginv(u.var.mod, eps=epsilon)

    dof.cov.model <- u.var.mod.ginv$rank
    test.cov.model <- (u.noint)%*%u.var.mod.ginv$Ginv%*%t(u.noint)

    # constrained
    test.cov.model.cons <- (u.noint.cons)%*%u.var.mod.ginv$Ginv%*%t(u.noint.cons)

    #### COVARIATE EFFECT TEST WITH ROBUST VARIANCE ####
    # not constrained
    u.var.robust <- sum.u2u2 - (sum.u2u1)%*%(1/sum.u1u1)%*%t(sum.u2u1)
    u.var.robust.ginv <- Ginv(u.var.robust, eps=epsilon)
    dof.cov.robust <- u.var.robust.ginv$rank
    test.cov.robust <- (u.noint)%*%u.var.robust.ginv$Ginv%*%t(u.noint)

    # constrained
    test.cov.robust.cons <- (u.noint.cons)%*%u.var.robust.ginv$Ginv%*%t(u.noint.cons)

    
    ## PREPARE RETURN OBJECT
    obj <- list(test.linkage=test.linkage,
              test.linkage.cons=test.linkage.cons,
              test.linkcov=test.linkcov,
              dof.linkcov=dof.linkcov,
              test.linkcov.cons=test.linkcov.cons,
              dof2.linkcov.cons=dof2.linkcov.cons,
              dof1.linkcov.cons=dof1.linkcov.cons,
              test.cov.model=test.cov.model,
              dof.cov.model=dof.cov.model,
              test.cov.model.cons=test.cov.model.cons,
              test.cov.robust=test.cov.robust,
              dof.cov.robust=dof.cov.robust,
              test.cov.robust.cons=test.cov.robust.cons,
              z.ped = z.ped)
  } else {  # if ncov>1
    obj <- list(test.linkage=test.linkage,
                test.linkage.cons=test.linkage.cons,
                z.ped = z.ped)
    
  }  #ncov=1

  return(obj)
}
