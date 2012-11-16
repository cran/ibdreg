#$Author: sinnwell $
#$Date: 2006/05/17 20:57:40 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.linkage.all.q,v 1.2 2006/05/17 20:57:40 sinnwell Exp $
#$Locker:  $
#$Log: print.linkage.all.q,v $
#Revision 1.2  2006/05/17 20:57:40  sinnwell
#remove unexpected sharing print
#
#Revision 1.1  2006/03/08 16:41:06  sinnwell
#Initial revision
#

######################################
## Jason Sinnwell
## Mayo Clinic, Div. of Biostatistics
## 11/2005
######################################

print.linkage.all <- function(x,
                              digits=max(options()$digits - 2, 5),
                              ...) {

  if(class(x)[1] != "linkage.all") stop("x must be a linkage.all object")
  test.df.names <-  c("pos(cM)", "Score Test", "d.f.", "pvalue")
  all.frame <- minpRows(obj=data.frame(x$all.frame[,1:2],
                          apply(x$all.frame[,3:6], 1, paste, collapse=":"),
                          x$all.frame[,7]),
                        colnames=test.df.names, rowname="constrained Linkage")

  
  printBanner(paste(x$status.method, "PAIRS", sep=' '))
  
  cat("Number of pedigrees used: \t", format(c(x$npeds, 123456))[1], "\n")
  cat("Number of relative pairs: \t", format(c(x$npairs, 123456))[1], "\n")
  cat("                AA pairs: \t", format(c(x$nstatus$AA, 123456))[1], "\n")
  cat("                UU pairs: \t", format(c(x$nstatus$UU, 123456))[1], "\n")
  cat("                AU pairs: \t", format(c(x$nstatus$AU, 123456))[1], "\n\n") 

  printBanner("Score Test for Linkage")

  print.data.frame(all.frame, digits=digits, ...)
  cat("\n\n")

  invisible(x)

}
