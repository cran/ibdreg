#$Author: sinnwell $
#$Date: 2006/05/19 18:14:32 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/plotpval.q,v 1.3 2006/05/19 18:14:32 sinnwell Exp $
#$Locker:  $
#$Log: plotpval.q,v $
#Revision 1.3  2006/05/19 18:14:32  sinnwell
#fix indexing to lwd
#
#Revision 1.2  2006/05/05 20:24:16  sinnwell
#add lwd, and change type to more options
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

plotpval <- function(pos, pmat, lty=1, lwd=1, col=1, title="", legend,
                     xlab="Position(cM)", ylab="-log10(pvalue)") {

  ## plot the -log10(pval) of pvalues, which are columns in pmat
  ## pos are chromosome positions
  ## add a legend with the same lines

  ncol.pmat <- ncol(as.matrix(pmat))

  ## if lty, lwd, or col are not same length as ncol(pmat)
  lty <- if(length(lty)==ncol.pmat) lty else rep(lty[1], ncol.pmat)
  col <- if(length(col)==ncol.pmat) col else rep(col[1], ncol.pmat)
  lwd <- if(length(lwd)==ncol.pmat) lwd else rep(lwd[1], ncol.pmat)
  
  ## find the range based on the minimum pvalue,  
  minp <- min(c(.001,pmat))  # .001 corresponds to upper ylim=3
  ylim=c(0,-log10(minp))
  many.pval <- ncol.pmat > 1

  #plot the first line
  y1 <- if(many.pval) -log10(pmat[,1]) else -log10(pmat)
  
  plot(x=pos, y=y1, ylim=ylim, type='l',
       xlab=xlab, ylab=ylab, lty=lty[1], col=col[1], lwd=lwd[1])

  title(title)

  # if more cols of pmat, plot them against pos
  if(many.pval) {
    for(k in 2:ncol.pmat) {
      y.new <- -log10(pmat[,k])
      lines(x=pos, y=y.new, lty=lty[k], lwd=lwd[k], col=col[k])
    }
    # add a legend
    legend(x=c(min(pos), min(pos)+.2*max(pos)), y=c(-log10(minp)*.75, -log10(minp)),
           legend=legend, 
           lty=lty, col=col, lwd=lwd, bty="n")
  }
  
}
