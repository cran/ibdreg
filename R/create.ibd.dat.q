#$Author: sinnwell $
#$Date: 2008/10/14 16:46:05 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/create.ibd.dat.q,v 1.11 2008/10/14 16:46:05 sinnwell Exp $
#$Locker:  $
#$Log: create.ibd.dat.q,v $
#Revision 1.11  2008/10/14 16:46:05  sinnwell
#change T to TRUE
#
#Revision 1.10  2006/10/30 16:49:43  sinnwell
#refine checks for post2 and prior2 ==0 for Xlinked
#
#Revision 1.8  2006/09/12 21:14:37  sinnwell
#merge.ibd.sex to mergeIBD for S3 check
#
#Revision 1.7  2006/07/06 19:05:40  sinnwell
#merge.ibd.sex function call change
#
#Revision 1.6  2006/07/06 14:51:07  sinnwell
#add x.linked options, which switch post and prior for male pairs from merlin
#
#Revision 1.5  2006/04/20 19:21:14  sinnwell
#postfile and priorfile parameter names
#
#Revision 1.4  2006/04/18 16:18:22  sinnwell
#null to prior
#
#Revision 1.3  2006/03/15 17:16:18  sinnwell
#a to prior
#
#Revision 1.2  2006/03/15 15:48:23  sinnwell
#change P1,P2 to post1,post2
#
#Revision 1.1  2006/03/08 16:46:01  sinnwell
#Initial revision
#
#####################################
## Jason Sinnwell Dan Schaid
## 11/2005
## Mayo, Division of Biostatistics
#####################################

create.ibd.dat <- function(postfile,          # ibd file with all loci, all alleles
                        priorfile,         # ibd file from one locus, no alleles
                        software="merlin", # ibd software, format differs 
                        x.linked=FALSE,     # if chrom 23, x.linked
                        cov.data=NULL,     # covariate data frame with sex, and person.id
                        rm.noninform=TRUE  # remove non-informative pairs?
                        ) {
 
  # check rm.noninform for keeping non-informative pairs
  if(missing(rm.noninform)) rm.noninform <- TRUE
  
  # combine a merlin output file which has FAMILY, ID1, ID2, post0,post1,post2
  # with merlin output run on a prior locus, those are prior ibd probs
  
  dat.ibd <- read.table(file=postfile, header=TRUE)

  # to make steps go faster, remove self rows, obviously not useful


 # SORT IS DONE WITHIN IBD.DF.MERLIN()
 # sort the data.frame by family, id1, id2, marker
 #  ord <- order(dat.ibd[,1], dat.ibd[,2], dat.ibd[,3])
 #  dat.ibd.sort <- dat.ibd[ord,]

 # construct post0, post1, post2 data.frames, columns are ibd prob by position
 # keep with ped, per1, per2 id's
  if(software=="merlin") {

    post.ibd.df <- ibd.df.merlin(dat.ibd)
    
    # here we modified ibd.df.merlin.q to accept merlin output column names
    # we change names to ped.id, person1.id, person2.id, post0, post1, post2

  } else { stop("software not recognized") }

  
 # Assume ibd software was run for Prior IBD, using 1 marker, no alleles
 # applies to use for all positions, all chromosomes
 # --exception for X has a different IBD file

  dat.prior <- read.table(file=priorfile, header=TRUE)
  dat.prior <- dat.prior[,-4]
  names(dat.prior) <- c("ped.id", "person1.id", "person2.id", "prior0", "prior1", "prior2")

  tmp1 <- dat.prior$person1.id
  tmp2 <- dat.prior$person2.id
  dat.prior$person1.id <- ifelse(tmp1 < tmp2, tmp1, tmp2)
  dat.prior$person2.id <- ifelse(tmp2 > tmp1, tmp2, tmp1)

  ord.prior <- order(dat.prior[,1], dat.prior[,2], dat.prior[,3])
  dat.prior <- dat.prior[ord.prior,]

 # merge Prior ibds with estimated ibds
  
  post.id.df <- data.frame(post.ibd.df[,1:3], indx.post=1:nrow(post.ibd.df))
  prior.id.df <- data.frame(dat.prior[,1:3], indx.prior=1:nrow(dat.prior))

  id.merge.df <- merge.data.frame(post.id.df, prior.id.df, by=1:3)

  ibd.df <- data.frame(post.ibd.df[id.merge.df$indx.post,],
                                        dat.prior[id.merge.df$indx.prior, -c(1:3)])


  # remove uninformative pairs
  epsilon <- .000001

  self <- abs(ibd.df$prior2 - 1.0) < epsilon
  ibd.df <- ibd.df[!self,]

  # remove parent=child pairs
  parentChild <- abs(ibd.df$prior1 - 1.0) < epsilon
  ibd.df <- ibd.df[!parentChild,]

  if(x.linked) { 
  # For x.linked, remove male pairs with prior0 = 1
    exclude <- abs(ibd.df$prior0 - 1.0) < epsilon
    ibd.df <- ibd.df[!exclude,]
  }


 # subset to informative affected relative pairs
  if(rm.noninform) {
    non.inform <- apply((ibd.df$post0 - ibd.df$prior0) < epsilon, 1, all) &
    apply((ibd.df$post1 - ibd.df$prior1) < epsilon, 1, all) &
    apply((ibd.df$post2 - ibd.df$prior2) < epsilon, 1, all)
    
    ibd.df <- ibd.df[!non.inform,]
  }
  
  if(x.linked && software=="merlin") {
    # add the sex covariate
    dat <- cov.data[c("ped.id", "person.id", "sex")]
    ibd.df <- mergeIBD(ibd.df, dat)

    male.pair <- (ibd.df$sex.1 == 1) & (ibd.df$sex.2 == 1)

    if(any(male.pair))
      {
        # interchange priors if male pairs

        prior2 <- ibd.df$prior2
        prior1 <- ibd.df$prior1

        change <- male.pair & (prior2 > epsilon) & (prior1 < epsilon)

        if(sum(change)) {
          ibd.df$prior1[change] <- prior2[change]
          ibd.df$prior2[change] <- prior1[change]
        }

        # interchange posteriors if male pairs
        post2 <- ibd.df$post2
        post1 <- ibd.df$post1

        ck2 <- apply(post2 > epsilon, 1, any)
        ck1 <- apply(post1 < epsilon, 1, all)
        change <- male.pair & ck1 & ck2
        if(sum(change)) {
          ibd.df$post1[change,] <- post2[change,]
          ibd.df$post2[change,] <- post1[change,]
        }

        if( any(ibd.df$prior2[male.pair] > epsilon)) {
        
          ibd.df.err <- ibd.df[ibd.df$prior2[male.pair] > epsilon,,drop=FALSE]
          pairs <- apply(ibd.df.err[,c(1:3)], 1, paste, collapse="-")
          stop(cat("Male pairs \n", pairs, " \nhave prior2 != 0\n"))

        }
        
        if( any(apply(ibd.df$post2[male.pair,] > epsilon, 1, any))) {
          badPost2 <- apply(ibd.df$post2[male.pair,] > epsilon, 1, any)
          ibd.df.err <- ibd.df[male.pair,][badPost2,]
          pairs <- apply(ibd.df.err[,c(1:3)], 1, paste, collapse="-")
          stop(cat("Male pairs \n", pairs, " \nhave post2 != 0\n"))

        }
        
      }  
    
  }
  
  
  class(ibd.df) <- c("ibd.dat" , "data.frame")
  
  ibd.df
  
}
