#$Author: sinnwell $
#$Date: 2007/01/03 19:01:37 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/minpRows.q,v 1.4 2007/01/03 19:01:37 sinnwell Exp $
#$Locker:  $
#$Log: minpRows.q,v $
#Revision 1.4  2007/01/03 19:01:37  sinnwell
#re-add if statement in paste(), only places number if multiple rows
#
#Revision 1.3  2007/01/03 15:24:32  sinnwell
#re-entered deleted line of code assigning return.row
#
#Revision 1.2  2007/01/02 21:32:59  sinnwell
#enforce a non-null row.name parameter in data.frame
#
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
  
  if(!nrow(obj) || !n.cols) {
    warning("object has zero dimension")
    return()
  }
   
  if(length(colnames) != n.cols) colnames=1:n.cols

  indx <- which(obj[,col.indx]==min(obj[,col.indx], na.rm=TRUE))

  rowlabels <- if(!length(rowname)) 1:length(indx) else paste(rowname, if(length(indx)>1) 1:length(indx))

  return.row <- data.frame(obj[indx,], row.names=rowlabels)
   
  dimnames(return.row)[[2]] <- colnames
  return.row
  
}
