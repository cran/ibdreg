#$Author: sinnwell $
#$Date: 2006/11/07 21:34:23 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/ibdreg.q,v 1.7 2006/11/07 21:34:23 sinnwell Exp $
#$Locker:  $
#$Log: ibdreg.q,v $
#Revision 1.7  2006/11/07 21:34:23  sinnwell
#no.dom to nodom
#
#Revision 1.6  2006/10/25 14:49:51  sinnwell
#remove library(Matrix)
#
#Revision 1.5  2006/05/26 16:17:32  sinnwell
#add epsilon for svd cutoffs
#
#Revision 1.4  2006/04/18 21:18:43  sinnwell
#move parameter order
#
#Revision 1.3  2006/03/15 21:35:42  sinnwell
#check for prior1 instead of a1
#
#Revision 1.2  2006/03/15 15:46:29  sinnwell
#change P1, P2 to post1, post2
#
#Revision 1.1  2006/03/08 16:44:11  sinnwell
#Initial revision
#

########################################
#   Jason Sinnwell,  Daniel Schaid
#   Div of Biostatistics
#   Mayo Clinic, HSR  2005
########################################

ibdreg <- function(formula,
                   status.method, # AA/AU/UU/ALL, or lowercase
                   c.scale='nodom', # either 'nodom' (no-dominance)
                                     # or 'minimax' (lodpal)
                   data,
                   status,
                   ped.id,
                   person.id,
                   ibd.dat, # ibd.dat object
                   ibd.var, # ibd.var object
                   subset,
                   weights,
                   na.action,
                   min.pairs=1, # min pairs for analysis
                   epsilon=1e-5,# singular value cutoff in Ginv
                   ...)

{

## Overview
## -----------------------------------------------------------
## Quasi-Likelihood Score Statistics in a regression framework
## to test genetic linkage on affected relative pairs with covariates.  
## Perform the following tests on AA, AU, or UU relative pairs
## a) Linkage
## b) Linkage with covariate effect
## c) Covariate Effect, adjusted for linkage
##   --Model Based Variance
##   --Robust Variance
## These tests are both unconstrained and constrained for excess
##    allele sharing in one direction
## To test linkage on all relative pairs at once, 
## d) Linkage while constraining ibd sharing in the direction
##      where excess allele sharing is expected under linkage
  

  Call <- match.call()

  
# CHECK INPUT PARAMETERS
#--------------------------------------------------------- 
  # check classes of ibd.var and ibd.dat objects
  if(!sr.class(ibd.var)[1]=='ibd.var')
    stop(paste(deparse(substitute(ibd.var)), " must be an object with class ibd.var"))

  if(!sr.class(ibd.dat)[1]=='ibd.dat')
    stop("ibd.dat must be an object with class ibd.dat")
    ibd.dat.names <-  c("ped.id", "person1.id", "person2.id",
                            "post0", "post1", "post2", "prior0", "prior1", "prior2")
  match.names <- match(ibd.dat.names, names(ibd.dat))
  if(any(is.na(match.names))) stop("ibd.dat object is missing named element ",
                                   ibd.dat.names[is.na(match.names)], "\n")

  # match status.method, if not matched stop()
  status.method = casefold(status.method, upper=TRUE)
  status.code <- match(status.method, c("AA", "UU", "AU", "ALL"))
  if(is.na(status.code))
    stop("Unrecognized status.method, choose from: AA, UU, AU, or ALL")

# CHECK FORMULA, it should be one of these formats
# -------------------------------------------------------
# either: ~pairs.fun1(cov1,cov1) [+ pairs.fun2(cov2) ]
# or:     ~1
# report error if a response is given,
# other errors will be caught in model.frame setup
  ftrms.attr <- attributes(terms(formula))
  if(!(ftrms.attr$response==0))
    stop("'formula' should not have a response,
          response is posterior IBD prob. in ibd.dat object")
  

# EXTRACT NECESSARY VECTORS FROM data
#-----------------------------------------------------
#  status, person.id, ped.id
#  plus any covariate used in ~formula
#  also, anything needed by model.frame parameters
  
  # extract covariates from formula, or leave as ~1
  vars <- all.vars(formula)
  if(length(vars)==0) {  # intercept only
    
    # Need just one column to make model.frame execute
    #   with other args (e.g. na.exclude, status,..)
    # In S, ~1 will not set up column of 1s
    # Set-up column of ped.id, it is required in data, we ignore it later
   
    temp.formula <- as.formula(paste("~", deparse(substitute(ped.id)), sep=''))
  } else {
    temp.formula <- as.formula(paste('~', paste(vars, collapse='+'), sep=''))
  }

  temp.call <- call('model.frame', formula=temp.formula)

  # Steps to assign args to model.frame call, Terry Therneau 7/2005
  # assign desired parameters from Call to temp.call
  k=length(temp.call)

    # because R doesn't automatically name model.frame column names
    # let k keep track of index of last element of temp.call,
    # and assign the name (only when using R)

  for(i in c('data', 'weights', 'subset', 'na.action', 'ped.id', 'person.id', 'status'))
    {
      if(!is.null(Call[[i]])) {
        temp.call[[i]] <- Call[[i]]
        k <- k+1
        if(is.R()) names(temp.call)[k]=i
      }

    }
  # pull items in temp.call, search in 'data', then the rest of search path
  # na.action, subset and other model.frame operations complete
  temp.m <- eval(temp.call, sys.parent())


# CREATE DATA FRAME (PAIRS.FRAME), LIKE THIS:
# ------------------------------------------
#  _______________________________________________________________
#  id.df         | y.mat | status.df  | x.mat (adjusted by c.scale) 
#  pid id.1 id.2 | y.ibd | AA  UU  AU |  int(1)   cov1.....covk
#  --------------|-------|------------|------------------------------
#   1   1    2   |  s(k) |c*1  0   0  |   c*(1)   c*(x+y)  c*(x-y)^2
#   1   1    9   |       | 0  c*1  0  |  	
#  ----------------------------------------------------------------

# pairs.frame is actually a list of data.frames with the same nrows,
# with the exception id.rm.df
# We can extract: id.df, y.mat, status.df, x.mat, id.rm.df at once

  pairs.frame <- create.pairs.frame(ibd.dat, temp.m, formula, c.scale)


# CHECK ADEQUACY OF SAMPLE SIZES
#------------------------------------------------------------
#  must be large enough for asymptotic chi-square to hold
#  use status.tbl (created below)
  status.tbl <- data.frame(AA.count=sum(pairs.frame$status.df$AA>0),
                           UU.count=sum(pairs.frame$status.df$UU>0),
                           AU.count=sum(pairs.frame$status.df$AU>0))
  
  
  ALL.linkage <- AA.linkage <- UU.linkage <- AU.linkage <- NULL
  
# SWITCH ON STATUS.METHOD FOR WHICH ANALYSES TO RUN
#-----------------------------------------------------------------
  switch(status.method,
         
         "AA" = {     # CONCORDANT AA
             if(status.tbl$AA.count < min.pairs) {
               stop(paste("Minimum number of pairs", min.pairs, "not met by",
                          status.method, "pairs: ", status.tbl$AA.count, sep=" "))
             }
             ibd.var.AA <- align.ibd.var(pairs.frame$id.df[pairs.frame$status.df$AA>0,],
                                               ibd.var, epsilon=epsilon)

             AA.linkage <- linkage.tests(pairs.frame$y.mat[pairs.frame$status.df$AA>0,],
                                         pairs.frame$x.mat[pairs.frame$status.df$AA>0,,drop=FALSE],
                                         ibd.var.AA, status.method=status.method, epsilon=epsilon)

         },

         "UU" = {  #  CONCORDANT UU
             if(status.tbl$UU.count < min.pairs) {
                 stop(paste("Minimum number of pairs", min.pairs, "not met by",
                            status.method, "pairs: ", status.tbl$UU.count, sep=" "))
             }

             ibd.var.UU <- align.ibd.var(pairs.frame$id.df[pairs.frame$status.df$UU>0,],
                                            ibd.var, epsilon=epsilon)

             UU.linkage <- linkage.tests(pairs.frame$y.mat[pairs.frame$status.df$UU>0,],
                                         pairs.frame$x.mat[pairs.frame$status.df$UU>0,,drop=FALSE],
                                         ibd.var.UU, status.method=status.method,
                                         epsilon=epsilon)

         },

         "AU" = {  # DISCORDANT AU
             if(status.tbl$AU.count < min.pairs) {
               stop(paste("Minimum number of pairs", min.pairs, "not met by",
                          status.method, "pairs: ", status.tbl$AU.count, sep=" "))
             }
           
             ibd.var.AU <- align.ibd.var(pairs.frame$id.df[pairs.frame$status.df$AU>0,],
                                            ibd.var, epsilon=epsilon)

             AU.linkage <- linkage.tests(pairs.frame$y.mat[pairs.frame$status.df$AU>0,],
                                         pairs.frame$x.mat[pairs.frame$status.df$AU>0,,drop=FALSE],
                                         ibd.var.AU, status.method=status.method,
                                         epsilon=epsilon)

         },
         
         "ALL" = { # ALL PAIRS, ALL TESTS
           # a) ALL pairs
             ibd.var.all <- align.ibd.var(pairs.frame$id.df,
                                             ibd.var, epsilon=epsilon)
             ALL.linkage <- linkage.all(pairs.frame$y.mat, pairs.frame$status.df,
                                             ibd.var.all, epsilon=epsilon)
           
           # b) concordant AA pairs
             ibd.var.AA <- align.ibd.var(pairs.frame$id.df[pairs.frame$status.df$AA>0,],
                                            ibd.var, epsilon=epsilon)
             AA.linkage <- linkage.tests(pairs.frame$y.mat[pairs.frame$status.df$AA>0,],
                                         pairs.frame$x.mat[pairs.frame$status.df$AA>0,,drop=FALSE],
                                         ibd.var.AA, status.method="AA",
                                         epsilon=epsilon)

           # c) concordant UU pairs
             ibd.var.UU <- align.ibd.var(pairs.frame$id.df[pairs.frame$status.df$UU>0,],
                                            ibd.var, epsilon=epsilon)
             UU.linkage <- linkage.tests(pairs.frame$y.mat[pairs.frame$status.df$UU>0,],
                                           pairs.frame$x.mat[pairs.frame$status.df$UU>0,,drop=FALSE],
                                           ibd.var.UU, status.method="UU", epsilon=epsilon)

           # d) discordant AU pairs
             ibd.var.AU <- align.ibd.var(pairs.frame$id.df[pairs.frame$status.df$AU>0,],
                                             ibd.var, epsilon=epsilon)
             AU.linkage <- linkage.tests(pairs.frame$y.mat[pairs.frame$status.df$AU>0,],
                                         pairs.frame$x.mat[pairs.frame$status.df$AU>0,,drop=FALSE],
                                         ibd.var.AU, status.method="AU", epsilon=epsilon)
           } )  # switch


#[7] MAKE IBDREG OBJECT
#----------------------------------------------------

  relpair.tbl <- table(pairs.frame$id.df$ped.id)

  ibdreg.obj <- list(Call=Call,
                     relpair.tbl=relpair.tbl,
                     status.tbl=status.tbl,
                     data.rm.df=pairs.frame$data.rm.df,
                     pairs.rm.df=pairs.frame$pairs.rm.df,
                     AA.linkage=AA.linkage,
                     UU.linkage=UU.linkage,
                     AU.linkage=AU.linkage,
                     ALL.linkage=ALL.linkage)

  sr.class(ibdreg.obj) <- c("ibdreg", "list")
  
  return(ibdreg.obj)
  
}



