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

.First.lib <- function(lib, pkg) {
   if(is.R()) library.dynam("ibdreg", pkg, lib)
   
   
## below:extra commands from arp.gee that move perl scripts into R
## may be needed for getting Allegro or Genehunter output here
# get perl script lines into an object
# assign that object as a value in the arp.gee frame
#   gh.pl <- readLines(paste(lib, pkg,
#                       "exec/parse.genehunter.ibd.pl", sep="/"), n=-1)
#   merlin.pl <- readLines(paste(lib, pkg,
#                       "exec/parse.merlin.pl", sep="/"), n=-1) 
#   pos <- (1:length(search()))[search() == "package:arp.gee"]
 
#   assign("parse.genehunter.ibd.pl", gh.pl , pos = pos)
#   assign("parse.merlin.pl", merlin.pl , pos = pos)
   
 }
