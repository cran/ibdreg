#!/usr/local/bin/perl

# Create map and dat files for MERLIN from linkage par file

$narg=@ARGV;
($narg == 3) || die "usage: merlin.pl <prefile> <parfile> <chrom number>\n";

$prefile  = @ARGV[0];
$parfile  = @ARGV[1];
$chrom    = @ARGV[2];


open (FILEIN,  $parfile ) || die "Input  file could not be opened\n";

# Process 1st line of par file

$_ = <FILEIN>;
chomp;
@Fld = split(' ', $_, 9999);

$nloci     = $Fld[0];
$risklocus = $Fld[1];
$sexlinked = $Fld[2];

if($chrom==23 & $sexlinked!=1){
  print "Error: Par file not for sexlinked\n";
  exit;
}

if($chrom<23 & $sexlinked!=0){
  print "Error: Par file for sexlinked\n";
  exit;
}


#$program   = $Fld[3];

# Process 2nd line of par file

$_ = <FILEIN>;
chomp;
@Fld = split(' ', $_, 9999);
$mutlocus = $Fld[0];
$mutMale = $Fld[1];
$mutFemale = $Fld[2];
$haplofreq = $Fld[3];

# Skip 3rd line of par file
$_ = <FILEIN>;

# Process remaining lines for different locus types

for($loc=1;$loc<=$nloci;$loc++){
   $_ = <FILEIN>;
   chomp;
   @Fld = split(' ', $_, 9999);
   $loctype = $Fld[0];

   push(@locatrib, $loctype);

   push(@locname, $Fld[3]);

     # Process affection locus
   if($loctype==1){
      $_ = <FILEIN>;
      $_ = <FILEIN>;
      chomp;
      @Fld = split(' ', $_, 9999);
      $nliab = $Fld[0];

      $nskip = $nliab;
      if($sexlinked == 1){
	$nskip = 2 * $nliab;
      }

      for($i=1;$i<=$nskip;$i++){
         $_ = <FILEIN>; chomp;
         @Fld = split(' ', $_, 9999);
       }
    }

   # Process numbered (marker) locus
   if($loctype==3){
     $_ = <FILEIN>;
   }
 }

  # Skip line for sex difference, interference
  $_ = <FILEIN>;


  # Process line with recombination values
  $_ = <FILEIN>;
  chomp;
  @Fld = split(' ', $_, 9999);

# Finished with par file - close file
close(FILEIN);


# Create Merlin map file

open(MAP, "> merlin.map" ) || die "merlin.map  file could not be opened\n";

   # create list of recombination values

   @recval = splice(@Fld, 1, ($nloci-2));
   unshift(@recval,0.0);  # place 0.0 at front of recval, for first position of map

   printf(MAP "CHROMOSOME MARKER LOCATION\n");
   $cumdist = 0.0;

   for($i=1;$i<$nloci;$i++){
     $theta = @recval[$i-1];
     $cumdist = $cumdist + 100*(log( (1 + 2*$theta)/(1- 2*$theta) )/4); # Kosambi Map
     # Haldane map: $cumdist = $cumdist + 100*( -log(1-2*$theta)/2);

     printf(MAP "%s %s %8.2f\n",$chrom,@locname[$i],$cumdist);
   }

close(MAP);


# Create Merlin dat file
open(DAT, "> merlin.dat") || die "merlin.dat  file could not be opened\n";

   for($i=0;$i<$nloci;$i++){
     if(@locatrib[$i]==1) {
        printf(DAT "A %s\n",@locname[$i]);
        if($nliab>1) { printf(DAT "C Liab\n"); }
      }

      if(@locatrib[$i]==3) {
        printf(DAT "M %s\n",@locname[$i]);
      }
   }

close(DAT);

if($chrom==23)
  {
    #check for errors:
    $cmd = "minx -d merlin.dat -p ".$prefile." -m merlin.map --error > merlinerror.out";
    system($cmd);

    # wipe out error genotypes
    $cmd = "pedwipe -d merlin.dat -p ".$prefile." > merlinpedwipe.out";
    system($cmd);

    # run MERLIN w/o errors:
    $cmd="minx -d wiped.dat -p wiped.ped -m merlin.map  --grid 1  --ibd  > merlinWiped.out";
    system($cmd);
  }
else
  {
    # check for errors:
    $cmd = "merlin -d merlin.dat -p ".$prefile." -m merlin.map --error > merlinerror.out";
    system($cmd);

    # wipe out error genotypes
    $cmd = "pedwipe -d merlin.dat -p ".$prefile." > merlinpedwipe.out";
    system($cmd);

    # run MERLIN w/o errors:
    $cmd="merlin -d wiped.dat -p wiped.ped -m merlin.map --grid 1  --ibd > merlinWiped.out";
    system($cmd);
}

system("rm -f merlin.dat");
system("rm -f merlin.map");
system("rm -f merlinerror.out");
system("rm -f merlin.err");
system("rm -f merlinpedwipe.out");
system("rm -f merlinWiped.out");
system("rm -f wiped.ped");
system("rm -f wiped.dat");

