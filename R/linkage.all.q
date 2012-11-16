#$Author: sinnwell $
#$Date: 2006/05/17 19:55:04 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.all.q,v 1.6 2006/05/17 19:55:04 sinnwell Exp $
#$Locker:  $
#$Log: linkage.all.q,v $
#Revision 1.6  2006/05/17 19:55:04  sinnwell
#take out unexpected stats
#
#Revision 1.5  2006/04/26 15:46:15  sinnwell
#*** empty log message ***
#
#Revision 1.4  2006/04/26 14:04:01  sinnwell
#pass epsilon to linkage.all.pedloop
#
#Revision 1.3  2006/04/26 14:02:26  sinnwell
#add epsilon
#
#Revision 1.2  2006/04/25 15:51:23  sinnwell
#change label of Vu.u to ginvVu.u and fix z.ped.UU and z.ped.AU to correct order
#
#Revision 1.1  2006/03/08 16:37:23  sinnwell
#Initial revision
#

######################################
#   Jason Sinnwell,  Daniel Schaid
#   Div of Biostatistics
#   Mayo Clinic, HSR  2005
######################################

linkage.all <- function(
     y.mat,         #-estimated ibd for rel. pairs
     x.adj,         #-not a model matrix; 3 cols for AA, UU, AU pairs,
                    # adjusted by c.scale
     ibdvar.lst,    #-vector of lists, each element has:
	            # sv.ginv: gen-inverse of the variance of ibd estimates
		    # rank: rank of the generalized inverse
		    # n: no. people in the pedigree
     epsilon=1e-5   # cutoff for singular values in Ginv()
)
{

# 0) OVERVIEW
#  Give test statistics for tests on relpairs for overall linkage
#    constrain the Betas with a contrast matrix Rmat
#    if one of the pair types, AA, UU, AU has no relpairs, these tests will not work
#    save additional vectors to test for unexpected ibd sharing

  verbose <- FALSE   # verbose=TRUE # for debugging
  positions <- as.numeric(dimnames(y.mat)[[2]])
  npairs <- nrow(y.mat)
  nstatus <- data.frame(AA=sum(x.adj[,1]>0),
                        UU=sum(x.adj[,2]>0),
                        AU=sum(x.adj[,3]>0), row.names="pairs")
  npeds <- length(ibdvar.lst)

# 1) SET UP VECTORS FOR TEST RESULTS
#    For all of the tests, store the test statistic, d.f, and pvalue
#    in separate vectors.  Preset their size to number of positions

  # linkage test vectors
  npos <- ncol(y.mat)
  test.all <- pval.all <- numeric(npos)
  
  # contrast matrix for betas
  Rmat <- matrix(c(-1,1,0,
                   0,-1,0,
                   0, 0,1), ncol=3, byrow=TRUE)

#3) ITERATE OVER npos FOR TESTS IN ALL POSITIONS

  for(pos in 1:npos) {
    if(verbose) cat("pos: ", pos, "\n")

    # 3a) Loop over peds, get U and Vu
    save.pedloop <- linkage.all.pedloop(y.mat[,pos], as.matrix(x.adj),
                                        ibdvar.lst)
    
    # 3b) Compute T3, which applies the constraints of Rmat, using U and Vu above
    u.var.ginv <- Ginv(save.pedloop$u.var, eps=epsilon)$Ginv
    # Vu.u is beta.tilde from manuscript
    ginvVu.u <- u.var.ginv %*% matrix(save.pedloop$u.score, ncol=1)
    ls.save <- lsConstrain.fit(x=ginvVu.u,
                                b=rep(0,3),
                                s=u.var.ginv,
                                a=Rmat,
                                iflag=rep(0,3))
    
    ## check ls.save for convergence, other checks
    ## if failed, stop with error message
    if(ls.save$ifault != 0) {
      switch(as.character(ls.save$ifault),
             "1"=stop("lsConstrain, itmax exceeded."),
             "3"=stop("lsConstrain, invalid constrain function for some row ASA'=0"),
             stop("Unknown error in lsConstrain")
             )
    }

    min.dist <- ls.save$min.dist
    test.all[pos] <- matrix(save.pedloop$u.score,ncol=3)%*%ginvVu.u - min.dist

    # pvalue is from chibar distribution, weights from sigma
    sigma <- Rmat%*%u.var.ginv%*%t(Rmat)
    wgts <- chibar3.w(sigma)
    pval.all[pos] <- 1 - pchibar(test.all[pos], df=0:3, wt=wgts)

  } # for pos in npos

  
# PREPARE OBJECT FOR RETURN

  obj <- list(status.method='ALL', npairs=npairs,
              nstatus=nstatus, npeds=npeds,
              all.frame=data.frame(positions, test.all,
                dof0=rep(0,npos), dof1=rep(1,npos), dof2=rep(2,npos),
                dof3=rep(3,npos), pval.all=pval.all)
              )

  class(obj) <- "linkage.all"

  return(obj)

}
