#$Author: sinnwell $
#$Date: 2006/05/05 20:25:12 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.linkage.tests.q,v 1.3 2006/05/05 20:25:12 sinnwell Exp $
#$Locker:  $
#$Log: plot.linkage.tests.q,v $
#Revision 1.3  2006/05/05 20:25:12  sinnwell
#add lwd and more types
#
#Revision 1.2  2006/05/02 16:31:42  sinnwell
#remove model-var test with covariates, commented
#
#Revision 1.1  2006/03/08 16:34:12  sinnwell
#Initial revision
#

##################################
# Jason Sinnwell
# Daniel Schaid
# Mayo Clinic, HSR, Biostatistics
# 2005
##################################

plot.linkage.tests <- function(x, ...) {

  if(!inherits(x, "linkage.tests"))
    stop("Not a linkage.tests object")

  xpos <- x$linkage.frame$positions

  linkage.only <- x$ncov==1

  # pull off title parm from ellipsis, eg. an extra label for chrom, etc
  dots <- as.list(substitute(list(...)))[-1]
  title.indx <- match('title', names(dots))
  title.user <- if(!is.na(title.indx)) dots[title.indx]$title else NULL
  
  if(linkage.only) {

    minp <- min(c(.001,  # corresponds to upper ylim=3
                  x$linkage.cons.frame$pval.linkage.cons))
    title.txt <- paste(x$status.method, "PAIRS: Linkage Test")

    plotpval(pos=xpos, pmat=x$linkage.cons.frame$pval.linkage.cons,
             lty=1, col=1,lwd=2,
             title=paste(title.user, "\n", x$status.method,
                         "PAIRS: Constrained Linkage Test", sep=" "))

  } else {
    plotpval(pos=xpos,
             pmat=cbind(x$linkage.cons.frame$pval.linkage.cons,
                x$linkcov.cons.frame$pval.linkcov.cons,
                #x$cov.model.cons.frame$pval.cov.model,
                x$cov.robust.cons.frame$pval.cov.robust),
             lty=c(1,3,6), col=c(1,2,4), lwd=c(2,3,2),      #lty=c(1,3,5,6), col=1:4,
             title=paste(title.user, "\n", x$status.method,
                         "PAIRS: Constrained Linkage Tests with Covariates"),
             legend=c("Linkage without Covariates",
             "Linkage with Covariates",
             #"Covariate Eff, model-var",
             "Covariate Eff, robust-var"))

  }
  
  invisible(x)

}
