#$Author: sinnwell $
#$Date: 2006/09/14 19:26:43 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/mergeIBD.q,v 1.1 2006/09/14 19:26:43 sinnwell Exp $
#$Locker:  $
#$Log: mergeIBD.q,v $
#Revision 1.1  2006/09/14 19:26:43  sinnwell
#Initial revision
#


########################################
#   Jason Sinnwell,  Daniel Schaid
#   Div of Biostatistics
#   Mayo Clinic, HSR  2006
########################################


mergeIBD <- function(ibd.dat, sex.dat){

# merge ibd  and covar data with sex based on ped.id, person1.id, person2.id in ibd.dat,
# and ped.id, person.id in sex.dat

   ibd.dat <- as.data.frame(ibd.dat)

   nm.ibd.dat <- names(ibd.dat)
   if(sum(nm.ibd.dat=="ped.id")!=1)  stop("Need one columnin ibd.dat named 'ped.id'")
   if(sum(nm.ibd.dat=="person1.id")!=1) stop("Need one columnin ibd.dat named 'person1.id'")
   if(sum(nm.ibd.dat=="person2.id")!=1) stop("Need one columnin ibd.dat named 'person2.id'")

   nm.cov.dat <- names(sex.dat)
   if(sum(nm.cov.dat=="ped.id")!=1) stop("Need one column in sex.dat named 'ped.id'")
   if(sum(nm.cov.dat=="person.id")!=1) stop("Need one column in sex.dat named 'person.id'")
 

   id.ibd <- data.frame(as.numeric(ibd.dat$ped.id),
                        as.numeric(ibd.dat$person1.id),
                        as.numeric(ibd.dat$person2.id))
   names(id.ibd) <- c("ped.id","person1.id","person2.id")
   id.ibd <- cbind(id.ibd,rownum=1:nrow(id.ibd))


   id.dat <- data.frame(as.numeric(sex.dat$ped.id),
                        as.numeric(sex.dat$person.id))
   names(id.dat) <- c("ped.id","person.id")
   id.dat <- cbind(id.dat, indx=(1:nrow(id.dat)))

   m1 <- merge(id.ibd,id.dat,by.x=c(1,2),by.y=c(1,2),all.x=T)
   m1 <- m1[order(m1$rownum),]

   m2 <- merge(id.ibd,id.dat,by.x=c(1,3),by.y=c(1,2),all.x=T)
   m2 <- m2[order(m2$rownum),]

   indx1 <- m1$indx
   indx2 <- m2$indx

   x1 <- sex.dat[indx1, -pmatch(c("ped.id","person.id"),names(sex.dat)), drop=FALSE]
   x2 <- sex.dat[indx2, -pmatch(c("ped.id","person.id"),names(sex.dat)), drop=FALSE]

   nm <- names(x1)

   name1 <- names(ibd.dat)
   name2 <- paste(nm,1,sep=".")
   name3 <- paste(nm,2,sep=".")

   df <- data.frame(ibd.dat,x1,x2)
   names(df) <- c(name1,name2,name3)

   return(df)
}
