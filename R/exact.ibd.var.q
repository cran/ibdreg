#$Author: sinnwell $
#$Date: 2006/10/12 14:02:44 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/exact.ibd.var.q,v 1.1 2006/10/12 14:02:44 sinnwell Exp $
#$Locker:  $
#$Log: exact.ibd.var.q,v $
#Revision 1.1  2006/10/12 14:02:44  sinnwell
#Initial revision
#
#Revision 1.1  2006/04/19 22:14:50  sinnwell
#Initial revision
#

## Jason Sinnwell
## Mayo Clinic, Division of Biostatistics
## for S/R package: ibdreg
## 4/2006

exact.ibd.var <- function(file){

## read in a temporary output file from exact.ibdvar.pl
## to be made into an ibd.var object for SPLUS/R

  var.lst <- list()
  list.index <- 1
  
  x <- scan(file,what="")
  x.length <- length(x)
  x.index <- 1

  while(x.index < x.length){
    ibdvar.i <- list()

    #read pedid
    ibdvar.i$ped.id <- as.numeric(x[x.index])
    x.index <- x.index+1

    # read npairs
    npairs <- as.numeric(x[x.index])
    x.index <- x.index+1

    # read person1.id
    ibdvar.i$person1.id = as.numeric(x[x.index:(x.index+npairs-1)])
    x.index <- x.index + npairs

    # read person2.id
    ibdvar.i$person2.id = as.numeric(x[x.index:(x.index+npairs-1)])
    x.index <- x.index + npairs 

    # read mean IBD sharing vector
    ibdvar.i$sm = as.numeric(x[x.index:(x.index+npairs-1)])
    x.index <- x.index + npairs 

    # add on rows of var-covar matrix
    sv <- matrix(0, ncol=npairs)
    for (k in 1:npairs) {
      sv <- rbind(sv, as.numeric(x[x.index:(x.index+npairs-1)]))
      x.index <- x.index + npairs
    }
    
    # drop first column of zeros
    ibdvar.i$sv <- sv[-1,,drop=FALSE]

    var.lst[[list.index]] <- ibdvar.i
    list.index <- list.index+1
  }

  sr.class(var.lst) <- c("ibd.var", "list")
  return(var.lst)
  
}

