#$Author: sinnwell $
#$Date: 2006/11/20 20:42:47 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.ibd.var.q,v 1.3 2006/11/20 20:42:47 sinnwell Exp $
#$Locker:  $
#$Log: sim.ibd.var.q,v $
#Revision 1.3  2006/11/20 20:42:47  sinnwell
#F to FALSE
#
#Revision 1.2  2006/10/30 19:17:41  sinnwell
#nsim to n.sim
#
#Revision 1.1  2006/10/10 21:21:07  sinnwell
#Initial revision
#
## authors: Dan Schaid and Dan Folie
## purpose: simulate alleles by gene dropping for an ibd.var object
## package: ibdreg
## date:    Oct 2006

sim.ibd.var <- function(pedfile, male.code=1, female.code=2,
                        x.linked=FALSE,  n.sim=1000, print.resources=FALSE){

   dat <- read.table(pedfile, header=FALSE, row.names=NULL)

   uped <- unique(dat[,1])
   nped <- length(uped)

   save <- vector("list", nped)

   for(i in 1:nped){


       dat.tmp <- dat[dat[,1]==uped[i],]
   
       dat.tmp <- dat.tmp[,-1] # remove ped id

       person <- dat.tmp[,1]
       father <- dat.tmp[,2]
       mother <- dat.tmp[,3]
       sex <- dat.tmp[,4]

       if(!all(sex==male.code | sex==female.code)){
         stop("Invalid code for sex")
       }

       male <- ifelse(sex==male.code, TRUE, FALSE)
   
       affection <- dat.tmp[,5]

       ped <- data.frame(person=person, father=father, mother=mother, sex=sex, male, affection = affection, 
                         row.names=NULL)
       
       ped <- sim.ibd.setup(ped, x.linked=x.linked)

       save[[i]]$ped.id <- uped[i]


      # create indx for persons that will be paired for ibd
      # configuration

      # old code for affecteds only:  person.indx <- (1:nrow(ped))[ped$affection==2]
      # now uses all pairs of subjects

       person.indx <- (1:nrow(ped))

       if(print.resources){
         cat("\nResources for ped ",i,"\n")
         resources({
           simped <- sim.mark.prop(ped, x.linked=x.linked, n.iter = n.sim)
           tmp <- ibd.moments(simped, person.indx, male, x.linked=x.linked)
         })
       } else {
         simped <- sim.mark.prop(ped, x.linked=x.linked, n.iter = n.sim)
         tmp <- ibd.moments(simped, person.indx, male, x.linked=x.linked)
       }
   

       save[[i]]$person1.id <- tmp$person1.id
       save[[i]]$person2.id <- tmp$person2.id
       save[[i]]$sm <- tmp$sm
       save[[i]]$sv <- tmp$sv
     }   

   sr.class(save) <- c("ibd.var", "list")
   return(save)
 }
