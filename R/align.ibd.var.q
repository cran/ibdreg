#$Author: sinnwell $
#$Date: 2006/04/26 14:14:27 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/align.ibd.var.q,v 1.4 2006/04/26 14:14:27 sinnwell Exp $
#$Locker:  $
#$Log: align.ibd.var.q,v $
#Revision 1.4  2006/04/26 14:14:27  sinnwell
#epsilon
#
#Revision 1.3  2006/04/25 15:52:34  sinnwell
#add ped.id to return list, easier for debugging
#
#Revision 1.2  2006/04/05 20:30:24  sinnwell
#add line to remove NAs from match of ibd.dat on ibd.var family
#
#Revision 1.1  2006/03/08 16:43:32  sinnwell
#Initial revision
#

######################################
#   Jason Sinnwell,  Daniel Schaid
#   Div of Biostatistics
#   Mayo Clinic, HSR  2005
######################################

align.ibd.var <- function(id.df, ibd.var, epsilon=1e-5) {
  
  # SAVE ELEMENTS OF IBD.VAR IN A VECTOR LIST  
  # THAT ARE NEEDED IN QLSCORE TEST STATISTICS

  # id.df: ped.id/person1.id/person2.id
  # ibd.var: list of lists for each pedigree
  #      elements: ped.id, person1.id, person2.id, sm, sv

  upeds <- unique(id.df$ped.id)
  save.pedvar <- vector("list", length(upeds))

  for(k in 1:length(ibd.var)) {
  # elements of ped are subset of all peds in ibd.var
    this.ped <- ibd.var[[k]]$ped.id

    # find where id.df ped.id matches k-th ped.id from ibd.var
    # this will be the place in save that ibd.var info is stored
    match.uped <- c(1:length(upeds))[this.ped==upeds]
    
    if(length(match.uped)) {
    # get data.frame with ped.id, per1, per2 for matched ped
      id.df.sub <- id.df[id.df$ped.id==this.ped,]
    # paste pairs of person.id's
      pair.paste <- paste(id.df.sub$person1.id,id.df.sub$person2.id,sep='_')
      varpair.paste <- paste(ibd.var[[k]]$person1.id, ibd.var[[k]]$person2.id, sep='_')
    # match pairs from ibd.ped.dat and those informative ones with pedvar entries
      pair.match <- match(pair.paste, varpair.paste)
      pair.match <- pair.match[!is.na(pair.match)]

    # save the ginv, n=dim, and rank of the ibd.var$sv matrix
      ginv.tmp <- Ginv(ibd.var[[k]]$sv[pair.match,pair.match], eps=epsilon)

      if(is.na(ginv.tmp$rank)) {
        msg <- paste("Generalized inverse does not exist for ibd in ped.id = ", this.ped, "\n")
        warning(msg)
      }
      
      save.pedvar[[match.uped]]$n <- length(pair.match)
      save.pedvar[[match.uped]]$sv.ginv <- ginv.tmp$Ginv
      save.pedvar[[match.uped]]$rank <- ginv.tmp$rank
      save.pedvar[[match.uped]]$ped.id <- this.ped

    } # if !is.na(match.uped)
  }   # for ped

  return(save.pedvar)

}
