#$Author: sinnwell $
#$Date: 2006/11/21 21:01:47 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/zzz.ibdreg.q,v 1.2 2006/11/21 21:01:47 sinnwell Exp $
#$Locker:  $
#$Log: zzz.ibdreg.q,v $
#Revision 1.2  2006/11/21 21:01:47  sinnwell
#remove library(Matrix)
#
#Revision 1.1  2006/03/08 17:02:02  sinnwell
#Initial revision
#
#
#

#.onLoad <- function(lib, pkg) {
#   library.dynam("ibdreg", pkg, lib)
#}

#.onAttach <- function(lib, pkg) {
#   library.dynam("ibdreg", pkg, lib)
#}


.onUnload <- function(libpath) {
  library.dynam.unload("ibdreg", libpath)
}

