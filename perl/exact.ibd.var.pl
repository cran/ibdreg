#!/usr/local/bin/perl

## Jason Sinnwell and Dan Schaid
## Mayo Clinic, Division of Biostatistics
## Splus/R package: ibdreg
## 7/2006


## use a pre file and chromosome number to make variance-covariance matrices for ibd
## sharing values between sets of relative pairs
##
## create kmx files (from merlin) which make joint probabilities of all ibd possibilities
## for relpairs in a pedigree, given a homozygous "dummy" marker.
## This makes it possible to construct the expected mean IBD sharing for relpairs, and also
## the expectation of the products, both are needed for variance-covariance calculations
##
##  for chromosome 23, merlin treats males as homozygous, so any relative pair
##  with a male has incorrect ibd sharing values.  Change any 2->1 and 4->2 in any
##  relative pair that contains a male. The gender status will be used from the
##  .pre input file


$narg=@ARGV;
($narg == 3) || die "usage: exact.ibd.var.pl <prefile> chrom [1|23] <outfile>\n";

$chrom = @ARGV[1];  # chromosome number
$varfile=@ARGV[2];  # name out output file for variance-covariance data
$verbose=0;         # change to 1 to print details

## make a pre file for merlin to make null Prior ibd probabilities
## this sub-fuction exists in merlin.prior.ibd.pl
## but is different b/c of maleids for ch23
## make maleid.txt to store ids in each pedigree that are male
## needed to re-code merlin.kmx files, as stated above.  maleid.txt is always made,
## whether chrom==23 or not, so always remove it

@prenames = makeNullPreFile(@ARGV[0]);

# make a parameter file for merlin to make null Prior ibd probabilities
$parname = makeNullParFile($chrom);

my $prefile = $prenames[0];  #"tempNull.pre";
my $midfile = $prenames[1];  # maleid.txt
my $parfile  = $parname;      #"tempNull.par";

open (FILEIN,  $parfile ) || die "Input file could not be opened\n";

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

# run merlin with options for kmx file output (--ibd --matrices)

if($chrom==23){

  $cmd = "minx -d merlin.dat -p ".$prefile."   -m merlin.map --ibd  --matrices > merlin.out";
  system($cmd);

}
else{
  $cmd = "merlin -d merlin.dat -p ".$prefile." -m merlin.map --ibd --matrices > merlin.out";
  system($cmd);

}

# remove temporary files made for and by merlin, 
## remaining files are merlin.kmx
system("rm -f merlin.dat");
system("rm -f merlin.map");
system("rm -f merlin.out");
system("rm -f ".$prefile);
system("rm -f ".$parfile);


## Make a file that can be read into S/R as an ibd.var object from
## Merlin results (options: [--ibd --matrices]  makes a file in '.kmx' format
            ## Family 2
            ## Pairs 101-102 1-102 1-101
            ## Position 0.0000
            ## 002 0.0625  :: Prob(pairs 1,2 share 0; pair 3 shares 2)=.0625
            ## 211 0.125   :: Prob(pair 1 shares 2; pairs 2,3 share 1)=.125
            ## ...

## For each family, make these structures:
##    mean of ibd sharing for each relative pair
##    covariance matrices of the posterior ibd prob between
##    pairs of relative pairs

## Print results to a file that can be read by exact.ibd.var() in Splus/R,
## which uses scan()
##    pedID
##    npairs
##    id1 (sorted by <id1,id2>)
##    id2
##    mean(IBD) for each pair
##    var(IBD), the covariance matrix between ibd statistics between
##    sets of relative pairs


open(KMXFILE, "merlin.kmx") || die "could not open merlin '.kmx' file \n";

if($chrom==23) {
  open(MALEIDFILE, "$midfile");
  @maleids = <MALEIDFILE>;
  $mnum=0;
}

# read in the first line, assumed to be Family
$line = <KMXFILE>;
chomp($line);
@spline = split(/\s+/, $line);
$pedid = $spline[1];

$nextpedid = 0;


# set up output filehandle, we will print to it after each pedigree is processed
open(VARFILE, ">$varfile") || die "could not open ibd.var output file \n";

while (<KMXFILE>) {

  my $starttime = time;

  ## vectors id1, id2 for storing relpair ids, 
  ## probvec joint prob of the ibd sharing values stored in rows of Xmat
  my @id1 = ();
  my @id2 = ();
  my @probvec = ();
  my @Xmat = ();

  # process pairs line
  $line=$_;
  chomp($line);
  @splitline=split(/\s+/, $line);

  ## start an  indicator array for which pairs have male, chrom==23-only
  @pairhasmale = ();

  foreach $p (1 .. $#splitline) {
    #index from 1 to length to skip word: "Pairs"

    @idpair = split('\-', $splitline[$p]);

    if($chrom==23) {

      ## keep an indicator vector for which pairs have a male,
      ## for recoding alleles shared
      $p0 = $idpair[0]; # new scalers necessary b/c [0] subset
      $p1 = $idpair[1]; # has functional meaning in regrexpr of grep

      push(@pairhasmale, (grep(/^$p0|\s$p0\s|$p0$/, $maleids[$mnum]) ? 1 : 0) |
	   (grep(/^$p1|\s$p1\s|$p1$/, $maleids[$mnum]) ? 1 : 0 ));
    }

    if(@idpair[0] < @idpair[1]) {

      push(@id1, $idpair[0]);
      push(@id2, $idpair[1]);

    }
    else {
      push(@id1, $idpair[1]);
      push(@id2, $idpair[0]);
    }
  } # foreach $p

  # skip over position line
  $line = <KMXFILE>;

  ## read lines of ibd joint probability between all pairs of pairs
  while(<KMXFILE>) {
    chomp;
    my @splitline=split(/\s+/, $_);
    if(@splitline[0] eq "Family") {
      $nextpedid = @splitline[1];

      if($chrom==23){
	$mnum++;
      }
      last;
    }
    else {
      push(@probvec, $splitline[1]);

      # break apart $shareline[0], put in xmat matrix
      my @sharevec = split(//,$splitline[0]);

      # for x-linked, change merlin sharing with male pairs
      if($chrom == 23) {
	
	## because merlin treats males as homozygous for X-chromosome,
	## for any relpair having a male, set 2->1 and 4->1
	## for alleles shared ibd in merlin.kmx

	for my $m (0..$#sharevec) {
	  $sharevec[$m] = $pairhasmale[$m] > 0 ? $sharevec[$m] > 0 ? 1 : 0 : $sharevec[$m];
	}
      }

      if($verbose) {
	print "@sharevec\n";
      }

      push(@Xmat, \@sharevec);
    }
  } # end while(<KMXFILE>)

  # process one pedigree's ibdvar
  # return: vector of npairs, references to: id1[sorted], id2[sorted], mean, var
  my ($npairs, $id1Ref, $id2Ref, $meanRef, $varRef) = pedIBDvar((\@id1,\@id2, \@Xmat, \@probvec));


  # de-reference id1, id2 mean vectors, and variance vector-o-vectors
  my @id1 = @$id1Ref;
  my @id2 = @$id2Ref;
  my @mean = @$meanRef;
  my @var = @$varRef;

  # print output to OUTFILE for the pedigree
  print VARFILE "$pedid\n";
  print VARFILE "$npairs\n";
  print VARFILE "@id1\n";
  print VARFILE "@id2\n";
  print VARFILE "@mean\n";

  for my $k (0 .. ($npairs-1)) {
    for my $l (0 .. ($npairs-1)) {
      print VARFILE "$var[$k][$l] ";
    }
    print VARFILE "\n";
  }

  my $endtime = time;
  my $totaltime = $endtime - $starttime;

  if($verbose) {
    print "ped: $pedid, npairs: $npairs, time=$totaltime \n";
  }
  $pedid = $nextpedid;

}

close(VARFILE);

if($chrom==23) {
  close(MALEIDFILE);
}

system('rm -f merlin.kmx');
system('rm -f '.$midfile);




##  DEFINE SUB FUNCTIONS

sub makeNullParFile{

  ## set up parameter file for merlin on a Null marker
  ## return par file name

  my($chrom) = @_;
  $parname = "tempNull.par";
  open(OUTFILE, ">$parname") || die "Output  file could not be opened\n";

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

  return $parname;

}


sub makeNullPreFile{
  ## make a Null pre file from the input file name,
  ## It is a pre file for only the first marker

  ## Take off ped per dad mom sex mk1
  ##   which is the first 6 columns of pre files

  ## store in a temporary pre file, "tempNull.pre"

  ## for later use on chrom23, print male ids within each pedigree to file, maleid.txt

  my($filein)  = @_;
  my(@line) = ();
  my($i) = 0;

  my(@maleid) = ();
  my($ped) = 0;
  my($newped) = 0;
  my($prefile) = "tempNull.pre";
  my($malefile) = "maleid.txt";
  open (FILEIN,  $filein ) || die "Input file could not be opened\n";
  open (FILEOUT, ">$prefile")  || die "Output  file could not be opened\n";

  open (IDFILE, ">$malefile");
  while(<FILEIN>){

    chomp;
    @line = split;
    $newped = $line[0];
    if($ped != $newped) {
      if($ped > 0) {
	## print male id code for each pedigree to temp file
	print(IDFILE "@maleid\n");
      }
      $ped=$newped;
      @maleid = ();
    }

    if($line[4]==1) {
      push(@maleid, $line[1]);
    }

    for($i=0; $i<6; $i++){
      printf(FILEOUT "%s ",$line[$i]);
    }
    printf(FILEOUT "1 1\n");

  }

  # print the male id lines for the last pedigree
  print(IDFILE "@maleid \n");

  close(FILEIN);
  close(FILEOUT);
  close(IDFILE);

  return ($prefile, $malefile);

}


sub pedIBDvar{

# use joint ibd sharing probabilities within a pedigree to make
# a covariance matrix for ibd values between pairs of relative pairs

  my ($subID1Ref, $subID2Ref, $subXmatRef, $subProbvecRef) = @_;

  my @subid1 = @$subID1Ref;
  my @subid2 = @$subID2Ref;
  my @subXmat = @$subXmatRef;
  my @subProbvec = @$subProbvecRef;

  my $npairs = @subid1; # number of pairs

  # make and index of the correct sorted ids
  my @unsorted = (0 .. ($npairs-1));
  # sort the unsorted vector of indices by id1, then by id2
  my @sorted = sort {$subid1[$a] <=> $subid1[$b] || $subid2[$a] <=> $subid2[$b]} @unsorted;


  # make mean ibd sharing array for each relpair
  my $nx = $#subXmat;
  my @submean = (0) x $npairs;
  for my $n (0 .. ($npairs-1)) {
    for my $m (0..$nx) {
      $submean[$n] += $subXmat[$m][$sorted[$n]] * $subProbvec[$m];
    }
    # print "mu: $submean[$n] \n";
  }

  # build empty matrices: covar matrix and E(ibd1*ibd2):(Eprod)
  # array of arrays, all $npairs long. $npairs x $npairs
  my @Eprod = ();
  my @covmat = ();
  for my $k (0 .. ($npairs-1)) {
    my @zerovec = (0) x $npairs;
    push(@Eprod, \@zerovec);
    push(@covmat, \@zerovec);
  }

  # fill Eprod and covmat
  for my $k (0 .. ($npairs-1)) {
    for my $l ($k .. ($npairs-1)) {
      for my $i (0 .. $nx) {
	$Eprod[$k][$l] += $subXmat[$i][$sorted[$k]] * $subXmat[$i][$sorted[$l]] * $subProbvec[$i];

      }

      $covmat[$k][$l] = $Eprod[$k][$l] - $submean[$k]*$submean[$l];
      $covmat[$l][$k] = $covmat[$k][$l];

    }
  }
  my @sortedid1 = @subid1[@sorted];
  my @sortedid2 = @subid2[@sorted];

  # return all scalers; for arrays, the scaler reference (address)
  return($npairs, \@sortedid1, \@sortedid2, \@submean, \@covmat);
}

