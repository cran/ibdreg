#$Author: sinnwell $
#$Date: 2006/10/10 21:20:05 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/ibd.moments.q,v 1.1 2006/10/10 21:20:05 sinnwell Exp $
#$Locker:  $
#$Log: ibd.moments.q,v $
#Revision 1.1  2006/10/10 21:20:05  sinnwell
#Initial revision
#
## authors: Dan Schaid and Dan Folie
## purpose: calculate moments (mean, variance) for simulated ibd vectors
##          from sim.mark.prop on each pedigree
## package: ibdreg
## date:    Oct 2006

ibd.moments <- function(simped, person.indx, male, x.linked=FALSE){

  tmp <- expand.grid(person.indx,person.indx)
  tmp <- tmp[tmp[,2] < tmp[,1], ]
  ivec <- tmp[,2]
  jvec <- tmp[,1]

  ibd <- NULL

  for(i in 1:length(ivec)){

      a1 <- simped$mark1[ivec[i],]
      a2 <- simped$mark2[ivec[i],]
      b1 <- simped$mark1[jvec[i],]
      b2 <- simped$mark2[jvec[i],]

      if(x.linked)
      {

          if(male[ivec[i]] & male[jvec[i]]) # male-male pair
          {
              tmp <- 1*(a1==b1)

          } else if((!male[ivec[i]]) & (!male[jvec[i]])) # female-female pair
          {
              tmp <- 1*(a1==b1) + 1*(a1==b2) + 1*(a2==b1) + 1*(a2==b2)

           } else if((male[ivec[i]]) & (!male[jvec[i]]))  # male-female pair
           {
              tmp <- 1*(a1==b1) + 1*(a1==b2)
           } else if ((!male[ivec[i]]) & (male[jvec[i]]))  # female-male pair
           {
              tmp <- 1*(a1==b1) + 1*(a2==b1)
            }

      } else
      {
             tmp <- 1*(a1==b1) + 1*(a1==b2) + 1*(a2==b1) + 1*(a2==b2)
     }

      ibd <- cbind(ibd,tmp)

  }

  sv <-  var(ibd)
  sm <-  apply(ibd, 2, mean)
  dimnames(sv) <- NULL
  names(sm) <- NULL

  return(list(person1.id=simped$person[ivec], person2.id=simped$person[jvec], sm=sm,sv=sv))

}
