#$Author: sinnwell $
#$Date: 2006/11/28 21:16:39 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/linkage.tests.q,v 1.5 2006/11/28 21:16:39 sinnwell Exp $
#$Locker:  $
#$Log: linkage.tests.q,v $
#Revision 1.5  2006/11/28 21:16:39  sinnwell
#make dimnames character, R 2.4 on linux didn't like numeric colnames
#
#Revision 1.4  2006/05/26 15:21:00  sinnwell
#make ped.zscores a data frame with ped ids as colnames and pos as rownames
#
#Revision 1.3  2006/05/17 20:22:00  sinnwell
#include z.ped for ped sharing departures
#
#Revision 1.2  2006/04/26 14:13:45  sinnwell
#epsilon to pedloop for Ginv
#
#Revision 1.1  2006/03/08 16:37:23  sinnwell
#Initial revision
#

######################################
#   Jason Sinnwell,  Daniel Schaid   
#   Div of Biostatistics
#   Mayo Clinic, HSR  2006
######################################

linkage.tests <- function(
     y.mat,         # estimated ibd for rel. pairs, col names are positions
     x.adj,         # x model matrix, adjusted by c-scale
     ibdvar.lst,    # each element has:
	            # sv.ginv: gen-inverse of the variance of ibd estimates
		    # rank: rank of the generalized inverse
		    # n: num. people in the pedigree--dim of var matrix
     status.method,  # AA/UU- concordant, AU-discordant
     epsilon=1e-5
)
{
 
#  OVERVIEW
#  Give test statistics for ql score tests on relpairs
#  allow either concordant or discordant
#  tests are:
#	linkage w/o covariate, unconstrained and constrained
#	linkage with covariate, unconstrained and constrained
#	covariate effect, model-based variance
#            --and constrained
#	covariate effect, robust variance
#            --and constrained
  
  verbose <- FALSE   # verbose=TRUE # for debugging
  npairs <- nrow(y.mat)
  npeds <- length(ibdvar.lst)
  
  
# SET UP VECTORS FOR TEST RESULTS
#   For all of the tests, store the test statistic, d.f, and pvalue
#   in separate vectors.  Preset their size to number of positions

  ncov <- ncol(x.adj)
  npos <- ncol(y.mat)
  positions <- as.numeric(dimnames(y.mat)[[2]])
  
  test.linkage <- test.linkage.cons <- numeric(npos)
  test.linkcov <- dof.linkcov <- numeric(npos)
  test.linkcov.cons <- dof2.linkcov.cons <- dof1.linkcov.cons <- numeric(npos)
  test.cov.model <- dof.cov.model <- numeric(npos)
  test.cov.model.cons <- test.cov.robust.cons <- numeric(npos)  
  test.cov.robust <- dof.cov.robust <- numeric(npos)
  test.cov.robust.cons <- dof.cov.robust.cons <- numeric(npos)
  z.ped.mat <- numeric()
  
## PRE-CALCULATE ALL BETA-0s
  # so just one loop needed over peds in making test stats
  B0.hat <- linkage.tests.B0(y.mat, xvec=x.adj[,1, drop=FALSE],
                       ibdvar.lst=ibdvar.lst)

#  ITERATE OVER npos FOR TEST STATS IN ALL POSITIONS
  for(pos in 1:npos) {
    if(verbose) cat("pos: ", pos, "\n")

    # Calculations are mostly done in a call a function that does pedloop, 
    #     which is in two forms for discordant and concordant
    #     The two functions differ only in constraint differences
    if(status.method=="AU") {
      save.pedloop <- linkage.tests.pedloop(y.mat[,pos], x.adj, ibdvar.lst,
                              B0=B0.hat[pos], concordant=FALSE, epsilon=epsilon)
    } else if(status.method=="AA" || status.method=="UU") {
      save.pedloop <- linkage.tests.pedloop(y.mat[,pos], x.adj, ibdvar.lst,
                              B0=B0.hat[pos], concordant=TRUE, epsilon=epsilon)
    } else {
      stop(paste("status method: ", status.method, " not recognized", sep=''))
    }
    
    # 3b) SAVE TEST STATS (AND DEG.FREEDOM) IN VECTORS

    # linkage only, unconstrained and constrained
    test.linkage[pos] <- save.pedloop$test.linkage
    test.linkage.cons[pos] <- save.pedloop$test.linkage.cons

    z.ped.mat <- rbind(z.ped.mat, save.pedloop$z.ped)
   
    
    if( ncov>1 ) { # only do these if covariate(s) given
      
      # test all betas==0, linkage with covariates
      dof.linkcov[pos] <- save.pedloop$dof.linkcov
      test.linkcov[pos] <- save.pedloop$test.linkcov
      # constrained linkcov, mixture chisq(q,q-1)
      dof2.linkcov.cons[pos] <- save.pedloop$dof2.linkcov.cons
      dof1.linkcov.cons[pos] <- save.pedloop$dof1.linkcov.cons
      test.linkcov.cons[pos] <- save.pedloop$test.linkcov.cons

      # test covariate effect with model-based variance
      dof.cov.model[pos] <- save.pedloop$dof.cov.model
      test.cov.model[pos] <- save.pedloop$test.cov.model
      # constrained cov.model
      test.cov.model.cons[pos] <- save.pedloop$test.cov.model.cons

      # test cov effect with robust variance
      dof.cov.robust[pos] <- save.pedloop$dof.cov.robust
      test.cov.robust[pos] <- save.pedloop$test.cov.robust
      # constrained cov.robust
      test.cov.robust.cons[pos] <- save.pedloop$test.cov.robust.cons

    } # if ncov>1

  } # for pos in npos

# GIVE DIMNAMES TO DATA FRAME WITH ZSCORES FOR PEDIGREES
# ROWS=POSITIONS, COLUMNS = PED.ID
  peds <- unlist(lapply(ibdvar.lst, function(lst) { lst$ped.id } ))
  z.ped.df <- as.data.frame(z.ped.mat, row.names=as.character(positions))
  names(z.ped.df) <- as.character(peds)
 
  
# CALCULATE P-VALUES FOR TEST STATISTICS
#    all non-constrianed tests are chi-square 
#    constrained with linkage are mixed chi-square (ncov, ncov-1),
#    constrained cov effect still chi-square(ncov-1)
#    where ncov = 1+rank(covar mat) for u scores of the pedigree

  dof.linkage <- rep(1,npos)
  dof1.linkage.cons <- rep(0,npos)
  dof2.linkage.cons <- rep(1,npos)

  pval.linkage <- 1-pchisq(test.linkage, dof.linkage)
  pval.linkage.cons <- 1-sapply(test.linkage.cons, pchibar,df=c(0,1), wt=c(.5,.5))
  
  if(ncov > 1) {

    pval.linkcov <- 1-pchisq(test.linkcov, dof.linkcov)

    # Mixture distribution deg. freedom for linkage & covar, constrained:
    #   dof2.linkcov.cons was set to dof.linkcov, test is mixture of dof2 and dof1=dof2-1.
    # If the unconstrained test had a cov matrix not full rank,
    #   then dof2 and dof1 would show it

    linkcov.cons.test.df <- cbind(test.linkcov.cons, dof1.linkcov.cons, dof2.linkcov.cons)
    pval.linkcov.cons <- apply(linkcov.cons.test.df,1,
                               function(vec)
                               1-pchibar(vec[1],df=vec[-1],wt=c(.5,.5)))

    pval.cov.model <-  1-pchisq(test.cov.model, dof.cov.model)

    pval.cov.model.cons <-  1-pchisq(test.cov.model.cons,dof.cov.model)

    pval.cov.robust <- 1-pchisq(test.cov.robust, dof.cov.robust)

    pval.cov.robust.cons <- 1-pchisq(test.cov.robust.cons, dof.cov.robust)

   
# PREPARE OBJECT FOR RETURN, DIFFERS FOR ncov=1

    obj <- list(status.method=status.method, B0.hat=B0.hat,
                npairs=npairs, npeds=npeds, ncov=ncov,
                linkage.frame = data.frame(positions, test.linkage,
                  dof.linkage, pval.linkage),
                linkage.cons.frame=data.frame(positions, test.linkage.cons,
                  dof1.linkage.cons, dof2.linkage.cons, pval.linkage.cons),
                linkcov.frame=data.frame(positions, test.linkcov,
                  dof.linkcov, pval.linkcov),
                linkcov.cons.frame=data.frame(positions, test.linkcov.cons,
                  dof1.linkcov.cons, dof2.linkcov.cons, pval.linkcov.cons),
                cov.model.frame=data.frame(positions, test.cov.model,
                  dof.cov.model, pval.cov.model),
                cov.model.cons.frame=data.frame(positions, test.cov.model.cons,
                  dof.cov.model, pval.cov.model.cons),
                cov.robust.frame=data.frame(positions, test.cov.robust,
                  dof.cov.robust, pval.cov.robust),
                cov.robust.cons.frame=data.frame(positions, test.cov.robust.cons,
                  dof.cov.robust, pval.cov.robust.cons),
                ped.zscore = z.ped.df
                )

    # end if ncov>1

  } else {  # ncov=1, returned object is different, only has linkage tests
    obj <- list(status.method=status.method, B0.hat=B0.hat,
                npairs=npairs, npeds=npeds, ncov=ncov,
                linkage.frame = data.frame(positions, test.linkage,
                  dof.linkage, pval.linkage),
                linkage.cons.frame=data.frame(positions, test.linkage.cons,
                  dof1.linkage.cons, dof2.linkage.cons,
                  pval.linkage.cons),
                ped.zscore = z.ped.df)
  }

  
  class(obj) <- c("linkage.tests", "list")

  return(obj)

}
