#$Author: sinnwell $
#$Date: 2006/11/21 17:26:24 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/print.linkage.tests.q,v 1.4 2006/11/21 17:26:24 sinnwell Exp $
#$Locker:  $
#$Log: print.linkage.tests.q,v $
#Revision 1.4  2006/11/21 17:26:24  sinnwell
#show.model.tests option enforced
#
#Revision 1.3  2006/11/09 21:47:04  sinnwell
#linkcov.cons.rows had linkcov.frame results, now fixed
#
#Revision 1.2  2006/03/31 21:35:08  sinnwell
#add show.mbtests option to print model-based. by default don't.
#make row names right-aligned for R. Splus already r-aligned.
#
#Revision 1.1  2006/03/08 16:41:06  sinnwell
#Initial revision
#

###################################
# Jason Sinnwell
# Daniel Schaid
# Mayo Clinic, HSR, Biostatistics
# 1/2006
###################################


print.linkage.tests <- function(x, digits = max(options()$digits - 2, 5), show.model.tests=FALSE, ...) {

  # For the tests that are part of linkage.tests object x, 
  # find index of position and test stats with smallest p-value for each test
  
  # Assemble rows of a return data.frame like this:
  #  "test-name" [position, test.stat, d.f., pval]
  #   where test-name is given as row name,
  #   and df is either numeric or the pasted columns if mixture d.f. 
  # Allow for multiple rows if two positions have the min pvalue

  # Tests include linkage-only
  # and additional tests with covariates if x$ncov > 1

  if(sr.class(x)[1] != "linkage.tests") stop("x must be a linkage.tests object")

  index <- 1:nrow(x$linkage.frame)

  test.df.names <- c("pos(cM)", "Score", "d.f.", "pvalue")
  bannertext <- "Score Tests for Linkage"

  linkage.only <- x$ncov==1
 
  # assemble linkage-only tests
  linkage.rows <- minpRows(x$linkage.frame, colnames=test.df.names, rowname="Linkage w/o Cov")

  linkage.cons.rows <- minpRows(obj=data.frame(x$linkage.cons.frame[,1:2],
            apply(x$linkage.cons.frame[,3:4], 1, paste, collapse=":"),
            x$linkage.cons.frame[,5]),
            colnames=test.df.names,rowname="constrained Linkage w/o Cov")

  test.df <- rbind.data.frame(linkage.rows, linkage.cons.rows)

  # now assemble other tests, if they exist
  if(!linkage.only) {
    bannertext <- paste(bannertext, " with Covariate(s)")

    linkcov.rows <- minpRows(x$linkcov.frame, colnames=test.df.names, rowname="Linkage w/ Cov")

    linkcov.cons.rows <- minpRows(obj=data.frame(x$linkcov.cons.frame[,1:2],
            apply(x$linkcov.cons.frame[,3:4], 1, paste, collapse=":"),
            x$linkcov.cons.frame[,5]),
            colnames=test.df.names,rowname="constrained Linkage w/ Cov")

    cov.model.rows <- minpRows(x$cov.model.frame, colnames=test.df.names,
                               rowname="Cov Effect (model)")

    cov.model.cons.rows <- minpRows(obj=x$cov.model.cons.frame,
                colnames=test.df.names, rowname="constrained Cov Effect (model)")

    cov.robust.rows <- minpRows(x$cov.robust.frame, colnames=test.df.names,
                               rowname="Cov Effect (robust)")

    cov.robust.cons.rows <- minpRows(obj=x$cov.robust.cons.frame,
                colnames=test.df.names, rowname="constrained Cov Effect (robust)")
    
    if(show.model.tests) {
      test.df <- rbind.data.frame(test.df,
                                  linkcov.rows,
                                  linkcov.cons.rows,
                                  cov.model.rows,
                                  cov.model.cons.rows,
                                  cov.robust.rows,
                                  cov.robust.cons.rows)
    } else {
      test.df <- rbind.data.frame(test.df,
                                  linkcov.rows,
                                  linkcov.cons.rows,
                                  cov.robust.rows,
                                  cov.robust.cons.rows)
      
    }
  }

  if(is.R()) {
    rnames <- rownames(test.df)
    maxwidth <- max(apply(matrix(rnames, ncol=1), 1, nchar))
    rownames(test.df) <- format(rnames, width=maxwidth, justify="right")
  }
   

  # print to the screen
  cat("\n")
  printBanner(paste(x$status.method, "PAIRS", sep=' '))
  maxdig=max(x$npeds, x$npairs)
  cat("Number of pedigrees used: \t", format(c(x$npeds, 123456))[1], "\n")
  cat("Number of relative pairs: \t", format(c(x$npairs,123456))[1],  "\n\n")

  printBanner(bannertext)
  print.data.frame(test.df, digits=digits, ...)

  cat("\n")
  
}


