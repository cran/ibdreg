#$Author: sinnwell $
#$Date: 2006/05/17 21:12:27 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plot.linkage.all.q,v 1.2 2006/05/17 21:12:27 sinnwell Exp $
#$Locker:  $
#$Log: plot.linkage.all.q,v $
#Revision 1.2  2006/05/17 21:12:27  sinnwell
#rm unexpect plot
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

plot.linkage.all <- function(x, ...) {

  if(!inherits(x, "linkage.all"))
    stop("Not a linkage.all object")

  positions <- x$all.frame$positions
  xpos <- as.numeric(positions)

  # extract user-given title from '...'
  dots <- as.list(substitute(list(...)))[-1]
  nm <- names(dots)
  title.indx <- match('title', nm)

  title.user <- if(!is.na(title.indx)) dots[title.indx]$title else NULL
  title.txt <- paste(title.user,
                "ALL PAIRS: Constrained Linkage Test ", sep='\n')

  minp <- min(c(.001,  # corresponds to upper ylim=3
                x$pval.all))
  ylim <- c(0,-log10(minp))

  linkage.y <- -log10(x$all.frame$pval.all)

  plot(x=xpos, y=linkage.y, ylim=ylim,
       xlab='Position', ylab='-log10(pvalue)',
       type='l')
  title(title.txt)
  
  invisible(x)

}
