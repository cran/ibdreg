#$Author: sinnwell $
#$Date: 2006/11/08 19:29:39 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/ibd.df.merlin.q,v 1.5 2006/11/08 19:29:39 sinnwell Exp $
#$Locker:  $
#$Log: ibd.df.merlin.q,v $
#Revision 1.5  2006/11/08 19:29:39  sinnwell
#undo model.matrix class, keep at matrix post0,1,2 frames
#
#Revision 1.4  2006/03/15 22:04:57  sinnwell
#change model.matrix names to post0, post1, post2
#
#Revision 1.3  2005/07/05 20:49:00  sinnwell
#expect names FAMILY, ID1, ID2, as from new merlin
#
#Revision 1.2  2005/01/13 20:38:24  sinnwell
#F to FALSE
#
#Revision 1.1  2004/10/07 21:51:38  sinnwell
#Initial revision
#
ibd.df.merlin <- function(ibd.dat){

  # Arguments: ibd.dat is a data.frame  The returned object is of class
  # ibd.df, which represents an "ibd data frame", where the ibd info
  # is stored in matrices, and these matrices are items in the returned
  # ibd.df


  # Create and return a dataframe of FAMILY, ID1, ID2, and ibd info
  # from merlin

  #  Change vec's to matrices of pairs (rows) by pos (cols)

  # count number of positions (for number of cols of post)
  zed <- ibd.dat$FAMILY==ibd.dat$FAMILY[1] & ibd.dat$ID1==ibd.dat$ID1[1] &
  ibd.dat$ID2==ibd.dat$ID2[1]
  
  npos <- sum(zed)

  npairs <- nrow(ibd.dat)/npos

  ord <- order(ibd.dat$FAMILY, ibd.dat$ID1, ibd.dat$ID2, ibd.dat$MARKER)
  ibd.dat <- ibd.dat[ord,]

  pos <- ibd.dat$MARKER[1:npos]

  ## model.matrix does not subset later in create.pairs.frame()
  ## It will work if keep these as matrix, then add to dat
  ## later like this: dat$post0<-post0

  post0 <- matrix(ibd.dat$P0,ncol=npos, byrow=TRUE)
  dimnames(post0) <- list(1:nrow(post0),pos)
#  oldClass(post0) <- c("model.matrix") #, "matrix")

  post1 <- matrix(ibd.dat$P1,ncol=npos, byrow=TRUE)
  dimnames(post1) <- list(1:nrow(post1),pos)
#  oldClass(post1) <- c("model.matrix") #, "matrix")

  post2 <- matrix(ibd.dat$P2,ncol=npos, byrow=TRUE)
  dimnames(post2) <- list(1:nrow(post2),pos)
#  oldClass(post2) <- c("model.matrix") #, "matrix")

  
  # Get ped,per1,per2 id's for first map position

  zed <- c(T,rep(F,(npos-1)))
  zed <- rep(zed, npairs)

  ped <-  ibd.dat$FAMILY[zed]
  tmp1 <- ibd.dat$ID1[zed]
  tmp2 <- ibd.dat$ID2[zed]
  per1 <- ifelse(tmp1 < tmp2, tmp1, tmp2)
  per2 <- ifelse(tmp2 > tmp1, tmp2, tmp1)

  # Calculate coefficients of ibd, and classify into pair types.
  # The order of pair type codes is from largest to smallest values of
  # P(ibd=1), with sibs higher than other types that have P(ibd=1)=0.5.
  # So, for example, sibs have type=1, avuncular type = 2, cousins type = 3, etc.

  # keep post0,1,2 as matrices in dat.  JPS for ibdreg, 11/2006
  # Splus couldn't subset when they were model.matrix class
  dat <- data.frame(ped.id=ped,person1.id=per1,person2.id=per2)
  dat$post0 <- post0
  dat$post1 <- post1
  dat$post2 <- post2

  
# data.frame methods need to work on ibd.df
 # sr.class(dat) <- c("ibd.df", "data.frame")

  return(dat)

}

