#$Author: sinnwell $
#$Date: 2006/05/17 20:58:13 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.ibdreg.q,v 1.2 2006/05/17 20:58:13 sinnwell Exp $
#$Locker:  $
#$Log: print.ibdreg.q,v $
#Revision 1.2  2006/05/17 20:58:13  sinnwell
#rm unexpected sharing print
#
#Revision 1.1  2006/03/08 16:41:06  sinnwell
#Initial revision
#


##################################
# Jason Sinnwell
# Daniel Schaid
# Mayo Clinic, HSR, Biostatistics
# 2005
##################################


print.ibdreg <- function(x,
                         digits = max(options()$digits - 2, 5),
                         ...) {
  cat("\n")
  print(x$Call)
  cat("\n")

  if(!is.null(x$AA.linkage))
    print.linkage.tests(x$AA.linkage, digits=digits, ...)

  if(!is.null(x$UU.linkage))
    print.linkage.tests(x$UU.linkage, digits=digits, ...)

  if(!is.null(x$AU.linkage))
    print.linkage.tests(x$AU.linkage, digits=digits, ...)

  if(!is.null(x$ALL.linkage))
    print.linkage.all(x$ALL.linkage,  digits=digits, ...)

  invisible(x)
}


