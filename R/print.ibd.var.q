#$Author: sinnwell $
#$Date: 2006/11/07 20:04:03 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.ibd.var.q,v 1.1 2006/11/07 20:04:03 sinnwell Exp $
#$Locker:  $
#$Log: print.ibd.var.q,v $
#Revision 1.1  2006/11/07 20:04:03  sinnwell
#Initial revision
#
#


##################################
# Jason Sinnwell
# Daniel Schaid
# Mayo Clinic, HSR, Biostatistics
# ibdreg package, 2006
##################################


print.ibd.var <- function(x,
                          ped.id = NULL,
                          sinkfile = NULL, 
                          digits = max(options()$digits - 2, 5),
                          ...) {

  # check ibd.var class
  if(!match("ibd.var", sr.class(x))) stop("Not an ibd.var object\n")

  
  ##-- if ped.id is NULL and sinkfile is NULL, don't print anything--TOO BIG
  ##-- if ped.id is not NULL, print the ibd.var elements that match ped.id's
  ##-- if sinkfile is not NULL, sink results to that file.
  
  if(!is.null(ped.id)) {
    
    if(!is.null(sinkfile)) sink(sinkfile)

    for(i in 1:(length(x))){

      if(!is.na(match(x[[i]]$ped.id, ped.id))) print(x[[i]], digits=digits)
      
    }

    if(!is.null(sinkfile)) sink()
    
  } else {  #ped.id not given, print all results to sinkfile
    
    if(is.null(sinkfile)) {
      xname <- deparse(substitute(x))
      sinkfile <- paste(xname, "sink", sep=".")
    }
    
    cat(paste("Sinking all pedigree ibd.var data to file", sinkfile, "...\n"))
    sr.class(x) <- "list"
    sink(sinkfile)
    print(x, digits=digits)
    sink()

  }  
  
  invisible(x)
  
}


