#$Author: sinnwell $
#$Date: 2006/11/28 21:14:31 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/create.pairs.frame.q,v 1.3 2006/11/28 21:14:31 sinnwell Exp $
#$Locker:  $
#$Log: create.pairs.frame.q,v $
#Revision 1.3  2006/11/28 21:14:31  sinnwell
#reorder ibd.df, y.mat, x.mat, statusid.df
#
#Revision 1.2  2006/03/15 16:21:06  sinnwell
#P=post, a=prior
#
#Revision 1.1  2006/03/08 16:46:25  sinnwell
#Initial revision
#
###########################################
## Jason Sinnwell and Dan Schaid
## 11/2005
## Mayo Clinic Division of Biostatistics
###########################################

create.pairs.frame <- function(ibd.dat,
                               model.mat,
                               formula,
                               c.scale)
{

  ## make three data.frames:
  ##     id.df [ ped.id, person1.id, person2.id]
  ##     y.mat [ estimated ibd sharing, rows: relpairs, cols: position]                
  ## status.df [ indicator columns for AA/UU/AU rel pairs, adjusted by c.vec] 
  ##     x.mat [ c.vec (intercept), cols for functions of covars for relpairs ]

  ped.id <- model.mat[['(ped.id)']]
  person.id <- model.mat[['(person.id)']]
  status <-  model.mat[['(status)']]

  # check values of status
  check <- 1*(status==1) + 1*(status==2)
  if(sum(check) != length(check))
    stop("Invalid value for status. Should be 1 (unaffected) or 2 (affected)")
  
  # model.mat has id and covariates for each person, rows removed for missing values
  # so match id1 and id2 from ibd.dat to those ids
  uid <- paste(ped.id, person.id, sep="_")
  id1 <- paste(ibd.dat$ped.id, ibd.dat$person1.id, sep="_")
  id2 <- paste(ibd.dat$ped.id, ibd.dat$person2.id, sep="_")
  match1 <- match(id1, uid, nomatch=NA)
  match2 <- match(id2, uid, nomatch=NA)

  # find persons in pedcov but not in ibd.dat
  # keep id, report error
  data.rm.id <- match(uid, c(id1,id2), nomatch=NA)

  if(any(is.na(data.rm.id))) {
    data.rm.df <- data.frame(ped.id=c(ibd.dat$ped.id, ibd.dat$ped.id)[is.na(data.rm.id)],
                 person.id=c(ibd.dat$person1.id, ibd.dat$person2.id)[is.na(data.rm.id)])
    data.rm.df <- data.rm.df[order(data.rm.df[,1]),]
#    cat("Persons in data not matched in ibd.dat object:\n")
#    print(data.rm.df)
#    stop("\n Please resolve \n")
  } else { data.rm.df <- NA }
  
  # find which pairs that have either person not matched in id,
  #   ie removed for NAs in the covariate matrix
  rm.pair <- apply(is.na(cbind(match1, match2)),1, any)
  
  # subset matched rows of id for pairs with none missing
  match1 <- match1[!rm.pair]
  match2 <- match2[!rm.pair]

  # MAKE DATA.FRAME WITH ped.id, person1.id person2.id
  id.df <- data.frame(ped.id=ibd.dat$ped.id[!rm.pair],
                      person1.id=ibd.dat$person1.id[!rm.pair],
                      person2.id=ibd.dat$person2.id[!rm.pair])

  # Make data.frame with ped.id person1.id person2.id of removed pairs
  pairs.rm.df <- data.frame(ped.id=ibd.dat$ped.id[rm.pair],
                      person1.id=ibd.dat$person1.id[rm.pair],
                      person2.id=ibd.dat$person2.id[rm.pair])
  
  ## MAKE C:SCALING FACTOR ON COVARIATES
  #  specific to ibd pair type, 
  c.vec <- create.pairs.frame.cvec(c.scale, prior1=ibd.dat$prior1[!rm.pair],
                              prior2=ibd.dat$prior2[!rm.pair])

  
  ## CREATE RESPONSE VECTORS, Y.MAT
  # the departure from expected ibd sharing, at each position
  y.mat <- (2*ibd.dat$post2[!rm.pair,]+ ibd.dat$post1[!rm.pair,]) -
           (2*ibd.dat$prior2[!rm.pair] + ibd.dat$prior1[!rm.pair])

  ## MAKE COLUMNS FOR COVARIATES
  # iterate over each term in call formula statement
  # create a function to evaluate, given the parameters subsetted for appropriate ids
  # evaluate and store column in x.mat

  x.mat <- matrix(c.vec, ncol=1)
  formula.terms <- attributes(terms(formula))$term.labels
  if(length(formula.terms)) {
    for(term in formula.terms) {

      ## split up the term by an open paren '(', use S::unpaste R::strsplit
      ## fixed=TRUE suggested 3/2013 by Claire Simpson
      split.term <- strsplit(term, split='(',fixed=TRUE)
      
      # extract the function name to apply to the covariate pair
      fname <- split.term[[1]][1]

      # either the parameter is double-listed or single
      # split by ',' or trim off ')'
      parname <- all.vars(formula(paste('~', term)))[1]

      # create a text string that is a call like: myfun(par[match1], par[match2])
      fcall=paste(fname, '(model.mat$', parname, '[match1], model.mat$', parname, '[match2])', sep='')

      # evaluate fcall to get xcol, center xcol
      xcol <-  eval(parse(text=fcall))
      xcol <- xcol-mean(xcol)
      xcol <- xcol*c.vec
      # add new column of x matrix the result in fcall, evaluated from m
      x.mat <- cbind(x.mat, xcol)

    } # for term

  } # if length(formula.terms)

  # set up 3 columns as indicators for AA, UU, and AU pairs
  # yet scaled down by c-scaling factor
  status1 <- status[match1]
  status2 <- status[match2]
  AA <- c.vec*((status1 + status2)==4)
  UU <- c.vec*((status1 + status2)==2)
  AU <- c.vec*((status1 + status2)==3)

  status.df <- data.frame(AA,UU,AU)
  
  # make sure id.df is ordered by pedigree, person1.id, person2.id
  # order other objects coinciding with it the same
  ord <- do.call("order", id.df)
  
  list(id.df=id.df[ord,], y.mat=y.mat[ord,, drop=FALSE], status.df=status.df[ord,], x.mat=x.mat[ord,, drop=FALSE],
       pairs.rm.df=pairs.rm.df, data.rm.df=data.rm.df)

}
