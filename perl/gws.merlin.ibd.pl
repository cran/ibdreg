#!/usr/local/bin/perl

## run a genome-wide scan using merlin
## which also creates prior ibd files for chrom 1 and 23
## and output files for creating ibd.var objects for chrom 1 and 23

  ## Assume that pre files are a directory named inputdir, and have the
  ## naming convention ch.k.pre for chromosome k=1..23

  ## Put posterior and prior ibd probability files in directories
  ## PostIBD and PriorIBD, respectively

## These assumptions may not be desired for users, so editing of 
## these directories and files may be necessary.

# Places that need editing have "Edit" comments above the neccessary lines

# consider a parameter to say whether to do exact ibd.var
# could also pass in the name of INPUTDIR as parameter
## $narg=@ARGV;
## ($narg == 3) || die "usage: gws.merlin.ibd.pl <indir> <outdir> <ibd.var option (1|0)>\n";

# Edit: define path to where pre and par files are stored (a single directory)
$dirInput = "INPUTDIR/";

# Edit: if exact ibd.var matrices desired set $doExactVar to 1
# the output file will be put in PriorIBD named "ibd.var.ch#.dat"
# where ch# is either 1 or 23
$doExactVar = 0;


# Edit: paths to where output files will be dumped
$dirPrior = "PriorIBD/";
$dirPost  = "PostIBD/";

system("mkdir ".$dirPrior);
system("mkdir ".$dirPost);

for($chrom = 1;  $chrom <= 23; $chrom++){

  # Edit: define file signature for pre file
  $prefile  = $dirInput."ch.".$chrom.".pre";

  # Edit: define file signature for par file
  $parfile  = $dirInput."ch.".$chrom.".par";

  # Note: use chrom=1 as template for all autosome prior IBDs and
  # prior kmx (merlin --matrices), and also run for chrom=23 for X
  # chromosome

  if(($chrom == 1) || ($chrom == 23)){

    system("merlin.prior.ibd.pl ".$prefile." ".$chrom);
    system("mv merlin.ibd ".$dirPrior."/chr".$chrom.".prior.ibd");

    ## create ibd.var.dat files to be made ibd.var objects by sr:exact.ibd.var()
    if($doExactVar) {
      system("exact.ibd.var.pl ".$prefile." ".$chrom." ibd.var.".$chrom.".dat");
      system("mv ibd.var."$chrom.".dat ".$dirPrior);
    }

  }

  # Posteriors for all chromosomes

  system("merlin.post.ibd.pl ".$prefile." ".$parfile." ".$chrom);

  system("mv merlin.ibd ".$dirPost."/chr".$chrom.".post.ibd");

}
