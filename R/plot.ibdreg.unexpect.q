#$Author: sinnwell $
#$Date: 2006/11/09 22:23:10 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.ibdreg.unexpect.q,v 1.1 2006/11/09 22:23:10 sinnwell Exp $
#$Locker:  $
#$Log: plot.ibdreg.unexpect.q,v $
#Revision 1.1  2006/11/09 22:23:10  sinnwell
#Initial revision
#
#Revision 1.1  2006/05/19 18:12:15  sinnwell
#Initial revision
#
##################################
# Jason Sinnwell
# Daniel Schaid
# Mayo Clinic, HSR, Biostatistics
# ibdreg package 2006
##################################

plot.ibdreg.unexpect <- function(x, status.method="AA", ...) {

  if(!inherits(x, "ibdreg"))
    stop("Not an ibdreg object")
  
  tests.obj <- switch(status.method,
                      "AA" = x$AA.linkage,
                      "AU" = x$AU.linkage,
                      "UU" = x$UU.linkage,
                      stop("status method not understood"))
  
  
  xpos <- tests.obj$linkage.frame$positions

  linkage.only <- tests.obj$ncov==1

  # pull off title parm from ellipsis, eg. an extra label for chrom, etc
  dots <- as.list(substitute(list(...)))[-1]
  title.indx <- match('title', names(dots))
  title.user <- if(!is.na(title.indx)) dots[title.indx]$title else NULL

  minp <- min(c(.001,  # corresponds to upper ylim=3
                tests.obj$linkage.cons.frame$pval.linkage.cons,
                tests.obj$linkage.frame$pval.linkage))

  plotpval(pos=xpos, pmat=cbind(tests.obj$linkage.cons.frame$pval.linkage.cons,
                       tests.obj$linkage.frame$pval.linkage),
           lty=c(1,2),  col=c(1,2), lwd=c(2,2),
           title=paste(title.user, "\n", tests.obj$status.method,
             "PAIRS: Constrained versus Unconstrained Linkage", sep=" "),
           legend=c("Constrained", "Unconstrained"))

  invisible(x)

}
