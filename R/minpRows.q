#$Author: sinnwell $
#$Date: 2006/03/08 16:42:40 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/minpRows.q,v 1.1 2006/03/08 16:42:40 sinnwell Exp $
#$Locker:  $
#$Log: minpRows.q,v $
#Revision 1.1  2006/03/08 16:42:40  sinnwell
#Initial revision
#

###################################
# Jason Sinnwell
# Daniel Schaid
# Mayo Clinic, HSR, Biostatistics
# 2/2006
###################################


minpRows <- function(obj, colnames=NULL, rowname=NULL, col.indx=ncol(obj)) {

  ## select a row from a data.frame for the minimum pvalue rows for test stats
  ##  -allow multiple rows to be selected
  ##  -allow row names to be given (repeat with numbers if multiple rows)
  ##  -allow assignment of column names

  n.cols <- ncol(obj)

  if(length(colnames) != n.cols) colnames=1:n.cols
  
  indx <- which(obj[,col.indx]==min(obj[,col.indx], na.rm=TRUE))

  rowlabels <- paste(rowname, if(length(indx)>1) 1:length(indx))

  return.row <- data.frame(obj[indx,], row.names=rowlabels)

  dimnames(return.row)[[2]] <- colnames
  return.row
  
}
