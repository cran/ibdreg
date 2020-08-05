#!/usr/local/bin/perl


$narg=@ARGV;
($narg == 2) || die "usage: merlin.ibd.null.pl <prefile> chrom [1-23]\n";


$chrom = @ARGV[1];
makeNullPreFile(@ARGV[0]);

makeNullParFile($chrom);

$prefile  = "tempNull.pre";
$parfile  = "tempNull.par";


open (FILEIN,  $parfile ) || die "Input  file could not be opened\n";

# Process 1st line of par file

$_ = <FILEIN>;
chomp;
@Fld = split(' ', $_, 9999);

$nloci     = $Fld[0];
$risklocus = $Fld[1];
$sexlinked = $Fld[2];
$program   = $Fld[3];

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
      if($sexlinked==1){
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



if($chrom==23){
  $cmd = "minx -d merlin.dat -p ".$prefile." -m merlin.map --ibd > merlin.out";
  system($cmd);
 # system("mv merlin.ibd merlinNullX.ibd");
}
else{
  $cmd = "merlin -d merlin.dat -p ".$prefile." -m merlin.map --ibd > merlin.out";
  system($cmd);
#  system("mv merlin.ibd merlinNull.ibd");
}

system("rm -f merlin.dat");
system("rm -f merlin.map");
system("rm -f merlin.out");
system("rm -f tempNull.pre");
system("rm -f tempNull.par");


sub makeNullParFile{
  my($chrom) = @_;

  open(OUTFILE, ">tempNull.par") || die "Output  file could not be opened\n";

  if($chrom==23)
  {
    printf(OUTFILE " 2 0 1 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM\n");
  }
   else
  {
    printf(OUTFILE " 2 0 0 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM\n");
  }

  printf(OUTFILE "0 0.0 0.0 0 << MUT LOCUS, MUT RATE, HAPLOTYPE FREQS (IF 1)\n");
  printf(OUTFILE "  1 2  << Order of loci\n");
  printf(OUTFILE "1   2  << disease locus\n");
  printf(OUTFILE " 0.997000 0.003000  << gene frequencies, wild type and mutant alleles\n");
  printf(OUTFILE " 1   << number of liability classes\n");
  printf(OUTFILE " 0.0010  1.0000 1.0000 << class 1\n");

  if($chrom==23){
    printf(OUTFILE " 0.0010  1.0000 << class  for males\n");
  }

  printf(OUTFILE "3 2 << Dummy-1\n");
  printf(OUTFILE "0.5 0.5 << dummy marker frequencies\n");
  printf(OUTFILE " 0 0 << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)\n");
  printf(OUTFILE " 0.5  << RECOMBINATION VALUES\n");
  printf(OUTFILE " 1 0.10000 0.45000     << REC VARIED, INCREMENT, FINISHING VALUE\n");

  close(OUTFILE);
}

sub makeNullPreFile{

  my($filein)  = @_;
  my(@line) = ();
  my($i) = 0;

  open (FILEIN,  $filein ) || die "Input  file could not be opened\n";
  open (FILEOUT, ">tempNull.pre")  || die "Output  file could not be opened\n";

  while(<FILEIN>){

    chomp;
    @line = split;

    for($i=0; $i<6; $i++){
      printf(FILEOUT "%s ",$line[$i]);
    }
    printf(FILEOUT "1 1\n");


  }

  close(FILEIN);
  close(FILEOUT);

  return;

}

sub isOdd{
  my($i) = @_;
  $odd = 1;
  if( (int($i/2)*2) == $i){
   $odd = 0;
 }
return($odd);
}
