#$Author: sinnwell $
#$Date: 2006/05/17 21:12:13 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.ibdreg.q,v 1.3 2006/05/17 21:12:13 sinnwell Exp $
#$Locker:  $
#$Log: plot.ibdreg.q,v $
#Revision 1.3  2006/05/17 21:12:13  sinnwell
#rm unexpect option
#
#Revision 1.2  2006/05/05 20:24:34  sinnwell
#give lwds for ALL and new type
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

plot.ibdreg <- function(x, ...) {

  # If status.method='ALL', and linkage-only
  #
  #      plot AA, UU, AU linkage lines on the same page
  #
  #      if unexpect=TRUE, another plot for those stats from ALL.linkage
  # else:
  #      plot of AA, UU, AU ALL linkage.tests obj, whichever is not null
  #         either a plot of one of AA, UU, AU, or
  #
  #         or all of them with ALL, and optional unexpect

  if(!is.null(x$ALL.linkage) && x$AA.linkage$ncov==1) {
    # plot AA, UU, AU, ALL lines on same plot
    dots <- as.list(substitute(list(...)))[-1]
    title.indx <- match('title', names(dots))
    title.user <- if(!is.na(title.indx)) dots[title.indx]$title else NULL

    plotpval(pos=as.numeric(x$AA.linkage$linkage.cons.frame$positions),
             pmat=cbind(x$ALL.linkage$all.frame$pval.all,
                        x$AA.linkage$linkage.cons.frame$pval.linkage.cons,
                        x$UU.linkage$linkage.cons.frame$pval.linkage.cons,
                        x$AU.linkage$linkage.cons.frame$pval.linkage.cons),
             lty=c(1,3,5,6), lwd=c(2,3,3,3), col=1:4,
             title=paste(title.user, "Constrained Linkage Tests by Status", sep="\n"),
             legend=c("ALL Pairs",
               "AA pairs",
               "UU pairs",
               "AU pairs"))

  } else {
    # plot whichever linkage results are not null
    if(!is.null(x$AA.linkage))
      plot.linkage.tests(x$AA.linkage, ...)
    if(!is.null(x$UU.linkage))
      plot.linkage.tests(x$UU.linkage, ...)
    if(!is.null(x$AU.linkage))
      plot.linkage.tests(x$AU.linkage, ...)
    if(!is.null(x$ALL.linkage))
      plot.linkage.all(x$ALL.linkage, ...)
  }

  invisible(x)
}
