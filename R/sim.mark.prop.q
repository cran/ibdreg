#$Author: sinnwell $
#$Date: 2006/11/21 19:58:22 $
#$Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.mark.prop.q,v 1.5 2006/11/21 19:58:22 sinnwell Exp $
#$Locker:  $
#$Log: sim.mark.prop.q,v $
#Revision 1.5  2006/11/21 19:58:22  sinnwell
#change F to FALSE on .C call
#
#Revision 1.4  2006/11/20 21:00:21  sinnwell
#add PACKAGE="ibdreg" to .C calls
#
#Revision 1.3  2006/11/07 19:50:19  sinnwell
#remove model.matrix caused problem in Splus, works in both S/R
#
#Revision 1.2  2005/10/04 18:42:48  folie
#Edited comments
#
#Revision 1.1  2005/09/28 18:37:57  folie
#Initial revision
#
# Function: sim.marker.prop()
# Author: Dan Folie
#
# This function takes an object "ped", of class "list" or "data.frame", that
# represents a pedigree data structure complete with markers for all of the
# founders present in the pedigree.  Output consists of an object of the same
# class as "ped" containing the simulated segregation of markers from the
# founders to the non-founder members of the pedigree.

sim.mark.prop <-
  function(ped,                # Simulated pedigree data - data.frame
           iseed = NULL,       # Random number generator seed 
           miss.val = c(NA,0), # Missing value codes 
           male.code = 1,      # Male sex code
           female.code = 2,    # Female sex code
           x.linked = FALSE,   # Tells if the locus is x.linked
           proband = 0,        # Starting proband (mainly used for testing)
           n.iter = 1)         # Number of times to propagate markers
{
  # Get random seeds for ruifAS183
  if(!is.null(iseed))
    { set.seed(iseed) }

  seed.array <- runif(3)

  # The seeds for the ranAS183 random number generator used in the C function
  # simulate_markers must be between 1 and 30000, but bigger is better
  # (we think), so we add 10000:

  iseed1 <- 10000 + 20000*seed.array[1]
  iseed2 <- 10000 + 20000*seed.array[2]
  iseed3 <- 10000 + 20000*seed.array[3]

  # Retrieve fields of ped that are to be used:
  person        <- ped$person
  father        <- ped$father
  mother        <- ped$mother
  sex           <- ped$sex
  chrom.markers <- cbind(ped$chrom1,ped$chrom2)


  # Proper sex not specified.
  if(any((sex != male.code) & (sex != female.code)))
    {
     if(!x.linked)
       warning("One or more subjects do not have there sex correctly specified.")
     else
       stop("One or more subjects do not have there sex correctly specified.")
   }
    
  # Check for males with heterozygous X-linked marker data:
  if(x.linked)
    {
     fail <- (chrom.markers[,1] != chrom.markers[,2]) & (sex==male.code)
     if(any(fail))
       {
        subjects <- person[fail]
        msg <- cat("The following males are heterozygous at an X-linked",
                   "marker:\n\n", subjects,"\n\n", sep=" ")
        stop(msg)
       }
    }
 
  if(is.factor(sex))
    {
     sex <- as.character(sex)
    }

  # Unspecified sex defaults to female
  sex[(sex != male.code) & (sex != female.code)] <- 0

  # Match male to sex code 1, and females to sex code 0
  sex[sex == male.code] <-   1
  sex[sex == female.code] <- 0

  # We will well-order the person, father, and mother vectors by considering
  # the entries of the person vector to be an ordering and then mapping this
  # ordering to 1,2,...,n where n = length(person):
  
  if(is.character(person))
    {     
     # All the codes involved in the pedigree should be present the in
     # person array:
     lev      <- unique(person)
     IDcodes  <- as.numeric(factor(c(father,mother,person),
                                   levels=lev,exclude=miss.val))
     IDcodes[is.na(IDcodes)] <- 0     
    }  
  
  if(is.factor(person))
    {
     # Convert factors to characters
     person  <- as.character(person)
     father  <- as.character(father)
     mother  <- as.character(mother)
     
     # All the codes involved in the pedigree should be present the in
     # person array:
     lev <- unique(person)
     IDcodes <- as.numeric(factor(c(father,mother,person),
                                  levels=lev,exclude=miss.val))
     IDcodes[is.na(IDcodes)] <- 0
    }  

  
  if(is.numeric(person))
    {     
     # All the codes involved in the pedigree should be present the in
     # person array:
     lev     <- unique(person)
     IDcodes <- as.numeric(factor(c(father,mother,person),
                                  levels=lev,exclude=miss.val))
     IDcodes[is.na(IDcodes)] <- 0
    }  
       
  n      <- length(IDcodes)/3
  
  father <- IDcodes[1:n]
  mother <- IDcodes[(n+1):(2*n)]
  person <- IDcodes[(2*n+1):(3*n)]
  nMark  <- ncol(chrom.markers)/2

  # Generate the pedigree:
  .C("sim_mark_prop_gen_ped",
     as.integer(person),
     as.integer(father),
     as.integer(mother),
     as.integer(sex),
     as.integer(n),
     as.integer(nMark),
     as.integer(chrom.markers),
     as.integer(iseed1),
     as.integer(iseed2),
     as.integer(iseed3),
     COPY = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
     PACKAGE="ibdreg")

  # Propagate markers n.iter times:
  mytemp <-
  .C("simulate_marker_propagation",
     as.integer(n.iter),
     ret=as.integer(integer(2*n.iter*n)),
     as.integer(proband),
     as.integer(x.linked),
     PACKAGE="ibdreg")$ret

  # Free up memory:
  .C("sim_mark_prop_free_mem",
     PACKAGE="ibdreg")

  # Construct return value:
  mytemp     <- matrix(mytemp, nrow=n)
  ## JPS changed 11/2/06 b/c of dimnames problem in S-PLUS
  ## keep mark1 and mar2 as matrix class
  #oldClass(mytemp) <- "model.matrix" 
  indices          <- 1:(2*n.iter)
  odd              <- as.logical(indices %% 2)
  even             <- !odd
  
  ped$mark1 <- mytemp[,odd, drop=FALSE]
  
  if(n.iter > 1)
     dimnames(ped$mark1) <- list(NULL, paste("m1.", 1:n.iter, sep =""))
  
  ped$mark2 <- mytemp[,even, drop=FALSE]
  if(n.iter > 1)
     dimnames(ped$mark2) <- list(NULL, paste("m2.", 1:n.iter, sep =""))

  ped$x.linked <- x.linked
  
  return(ped)
     
}
