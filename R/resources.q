## $Id: resources.q,v 1.1.1.1 2011/02/24 14:16:27 sinnwell Exp $
##

## S: Beth Atkinson and SR: Jason Sinnwell
## Mayo Clinic, Division of Biostatistics
## Under Strategic grant to make srlocal

# print resource details for the evaluation of an expression
# Reference, page 151-152 of S PROGRAMMING by Venebles and Ripley.

# example:
#resources({
#  norm <- rnorm(10000)
#  mat <- matrix(norm, ncol=100)
#  mat <- 2*mat
#})

resources <- function(expr, doPrint=TRUE) {

# Returns and prints the following information for the current S-PLUS session.
# Note: All times are expressed in units of seconds.
#   usertime - User elapsed times in S-PLUS.  
#   systemtime - System elapsed times in S-PLUS. 
#   cpu - The sum of usertime and systemtime.
#   elapsed - Elapsed time.
#   % cpu - (100 * cpu) / elapsed.
#   child - The sum of chusertime and chsystemtime.
#   usertimech - User elapsed times in child processes.
#   systemtimech - System elapsed times in child processes.
#   mem - maximum memory allocated by the process



  #loc <- sys.parent(1)
  #if(loc == 1) loc <- FALSE
  on.exit(cat("Timing stopped at:", proc.time() - time, "\n"))
  expr <- substitute(expr)
  stime <- proc.time()
  
  mem1 <- gc(reset=TRUE)
  w <- eval(expr)#, envir=loc) 

  etime <- proc.time()

  mem2=gc(reset=TRUE)
  mem <- data.frame(mem2[2,1]-mem1[2,1])
  names(mem) <- "Vcells"

  on.exit()
  
  if(length(stime)==3) stime <- c(stime,0,0)
  if(length(etime)==3) etime <- c(etime,0,0)
  time <- etime-stime
  time[3] <- max(time[3], time[1] + time[2])
  df <- data.frame(CPU = time[1] + time[2],
          Elapsed=time[3],
          CPU.pct=round((100*(time[1]+time[2]))/time[3],4),
          Child=time[4]+time[5],
          mem)
  if(doPrint) {
    print(df)
  }
  invisible(df)
    
}


