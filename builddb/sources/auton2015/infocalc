#!/usr/bin/perl
# Noah Rosenberg
# rnoah@umich.edu
# Infocalc - a program for calculating marker informativeness statistics
# Version 1.1 - December 23, 2006
# Version 1.0 - June 15, 2004
#
# New features in version 1.1
#   -Computation of the diploid version of ORCA from Rosenberg (2005)
#   -Incorporation of prior population weights on assignment probabilities
#   -More extensive command-line arguments
#
#########################################################################
#
# INSTRUCTIONS
# 
# This program computes the statistics I_n, I_a, and ORCA, as described in 
#           N A Rosenberg, L M Li, R Ward, J K Pritchard (2003) 
#           Informativeness of genetic markers for inference of 
#           ancestry.  Am J Hum Genet 73:1402-1422.
#
# The version of ORCA studied in Rosenberg et al. (2003) is based on 
# assignment of a single allele.  This program also computes a version of ORCA 
# based on assignment of a diploid genotype, as described in more detail in
#           N A Rosenberg (2005)
#           Algorithms for selecting informative marker panels for
#           population assignment.  J Comput Biol 12:1183-1201.
#
# When using infocalc, please cite the appropriate paper for the 
# statistic you are using.
#
# Each statistic is computed for each locus of a data file.  The data
# file must be in the 2-line format used by the program STRUCTURE, with
# a header line that lists the loci.  The main properties of this file
# format are:
#
# 1. The top line gives the list of loci.  The ith locus listed
#    corresponds to the data in the (i+5)th column.
# 2. On subsequent lines, the first five columns give individual diploid
#    genotypes.  In the example below, the five columns are individual
#    code, population code, population name, country, and geographic region.
# 3. Each individual is spread over two lines.  The two lines do not
#    indicate phase, and it is arbitrary whether an allele is placed on
#    the first or the second of the two lines.
# 4. Missing data is marked by a particular value, in this case -9.
# 
# An example data file, named mksp.stru, is below:
# 
# D9S1779 D9S1825 D7S2477 D17S784 D16S403 D3S1262 D10S189
# 854 86 Maya Mexico AMERICA 124 129 152 -9 138 112 186
# 854 86 Maya Mexico AMERICA 142 135 156 -9 140 124 186
# 855 86 Maya Mexico AMERICA 124 129 156 230 138 112 186
# 855 86 Maya Mexico AMERICA 124 129 164 234 140 112 186
# 856 86 Maya Mexico AMERICA 124 129 142 234 138 112 186
# 856 86 Maya Mexico AMERICA 124 135 156 234 146 114 188
# 995 82 Karitiana Brazil AMERICA -9 129 142 228 142 118 182
# 995 82 Karitiana Brazil AMERICA -9 129 156 234 148 124 182
# 996 82 Karitiana Brazil AMERICA 124 129 142 234 142 112 182
# 996 82 Karitiana Brazil AMERICA 124 137 156 234 142 124 182
# 997 82 Karitiana Brazil AMERICA 124 129 142 228 134 124 182
# 997 82 Karitiana Brazil AMERICA 146 135 142 228 144 124 182
# 998 82 Karitiana Brazil AMERICA 124 129 142 228 134 112 182
# 998 82 Karitiana Brazil AMERICA 146 135 142 234 142 124 182
# 999 82 Karitiana Brazil AMERICA 124 129 142 228 142 116 182
# 999 82 Karitiana Brazil AMERICA 146 129 156 228 144 124 182
# 830 83 Surui Brazil AMERICA 124 129 142 232 138 112 186
# 830 83 Surui Brazil AMERICA 124 135 156 234 146 124 186
# 832 83 Surui Brazil AMERICA 124 129 142 228 146 116 186
# 832 83 Surui Brazil AMERICA 124 129 156 232 146 124 186
# 833 83 Surui Brazil AMERICA 124 129 142 234 140 112 186
# 833 83 Surui Brazil AMERICA 144 135 156 234 146 112 186
# 834 83 Surui Brazil AMERICA 124 129 156 228 140 112 186
# 834 83 Surui Brazil AMERICA 144 129 156 228 146 112 186
# 835 83 Surui Brazil AMERICA 124 129 156 228 140 112 186
# 835 83 Surui Brazil AMERICA 124 135 156 228 150 112 186
# 837 83 Surui Brazil AMERICA 124 129 142 228 140 124 186
# 837 83 Surui Brazil AMERICA 124 129 156 234 146 124 186
# 838 83 Surui Brazil AMERICA 124 129 142 234 140 112 180
# 838 83 Surui Brazil AMERICA 124 135 156 234 146 112 186
# 1041 87 Pima Mexico AMERICA 124 129 142 234 140 112 180
# 1041 87 Pima Mexico AMERICA 124 129 142 234 150 112 182
# 1042 87 Pima Mexico AMERICA 124 129 142 230 150 112 180
# 1042 87 Pima Mexico AMERICA 152 129 142 234 150 124 182
# 1043 87 Pima Mexico AMERICA 124 129 142 232 142 112 180
# 1043 87 Pima Mexico AMERICA 146 135 142 234 142 112 182
# 1044 87 Pima Mexico AMERICA 124 129 142 232 150 112 182
# 1044 87 Pima Mexico AMERICA 124 129 142 234 150 112 182
#
# 
# RUNNING THE PROGRAM
# 
# The program is executed from a command line.  To run the program, type
#
#     ./infocalc
#
# This will run the program with default parameters, which can be
# modified with additional command-line arguments.
#
#     -column      Column to use as the population identifier (default 3)
#
#     -input       Input data file (default "infile")
#
#     -numpops     Number of populations in the data file (default 2)
#
#     -output      Output file for the results (default "outfile")
#
#     -weightfile  A file of population weights (see below - default is 
#                                                not to use a weightfile)
#
# The command-line options can be placed in any order. Also, the whole 
# flag does not need to be written, and the first  letter will suffice.  
# For example, the following two commands both will use column 3, input 
# file mksp.stru output file  mksp.stru.infostats, and weightfile 
# mksp.weight.
#
# ./infocalc -column 3 -numpops 4 -input mksp.stru -output 
#            mksp.stru.info-weighted -weightfile mksp.weight
#
# ./infocalc -n 4 -c 3 -o mksp.stru.info-weighted -i mksp.stru -w mksp.weight 
#
# The output is printed to mksp.stru.info-weighted.  From left to right, 
# the columns of the output file are the I_n, I_a, ORCA (1-allele),
# and ORCA (2-allele) statistics, as indicated in the header line.  
# There is one line of output for each locus.  
#
# In the output file, two lines follow the informativeness statistics 
# for the loci.  The first of these lines recapitulates the command used 
# to produce the output (with the command-line arguments placed in a 
# standardized order).  The second line indicates the weights used 
# for the populations, with each weight placed to the right of the
# corresponding population name.
#
# As described in the Rosenberg et al. (2003) paper, difficulties 
# sometimes arise in the calculation of I_a, so that denominators of 
# zero are produced.  In these cases, -9999  (or a number very close 
# to -9999, since some terms may be added to the -9999) is printed for 
# the value of I_a.  Also, when a weightfile is used or when more than
# 24 populations are specified, I_a is not calculated and "NA" is 
# printed in the I_a column.
#
# 
# THE INFORMATIVENESS STATISTICS, AND
# 1-ALLELE and 2-ALLELE VERSIONS OF ORCA 
#
# I_n (Rosenberg et al. 2003, eq. 4) gives the amount of information
# gained about population assignment from observation of a single
# randomly chosen allele at a locus.  This formula assumes a uniform
# prior across K potential source populations for the origin of the
# allele.
#
# I_a (Rosenberg et al. 2003, eq. 14) gives the amount of information
# gained about ancestry coefficients from observation of a single
# randomly chosen allele at a locus.  This formula assumes a uniform 
# prior over a simplex for the distribution of ancestry coefficient
# vectors in the population.
# 
# Eq. 10 of Rosenberg et al. (2003) gives the 1-allele version of ORCA.
# This is the optimal rate of correct assignment for a locus if a
# randomly chosen allele in a pooled collection of populations is
# assigned to its most likely source population.
# 
# A 2-allele version of ORCA is also possible.  This is the optimal
# rate of correct assignment for a locus if a randomly chosen diploid
# genotype is assigned to its most likely source population.  The
# 2-allele version of ORCA is obtained from the L=1 case of eq. 12 of
# Rosenberg et al. (2003), or from the M=1 case of eq. 1 of Rosenberg
# (2005).
#
# 
# WEIGHTS
# 
# The I_n and ORCA statistics are flexible in that they can accommodate
# a non-uniform prior on the source population of an individual (the I_a
# statistic does not have this flexibility).  This is accomplished by
# replacing 1/K in the appropriate equation with q_i, the prior
# probability that an individual is from population i (I_n - Rosenberg
# et al. 2003 eq. 4; ORCA - Rosenberg et al. 2003 eq. 10; 2-allele ORCA
# - Rosenberg et al. 2003 eq. 12 or Rosenberg 2005 eq. 1).
#
# Here is an example of a weightfile, mksp.weight.
#
# Maya 0.35
# Karitiana 0.15
# Surui 0.2
# Pima 0.3
#
# Running the command
#  ./infocalc -n 4 -c 3 -o mksp.stru.info-weighted -i mksp.stru -w mksp.weight 
# will measure the information content of markers subject to the prior
# that individuals have 0.35 probability of being from the Maya
# population, 0.15 from Karitiana, 0.2 from Surui, and 0.3 from Pima.
#
# Performing the analysis without a weightfile, as in the following,
#  ./infocalc -n 4 -c 3 -o mksp.stru.info -i mksp.stru
# will assume that each population is equally likely to be the source
# population.  The default is to assume that each source is equally likely.
# 
# The weightfile must contain the same population names as the data
# file, in any order.  The first column of the weightfile gives the
# names of the populations, and the second column gives the weights.
# The weights must sum to 1.  Additional whitespace between the columns
# is tolerated.  However, each row must begin with the population name
# (that is, no whitespace is allowed before the population name).
#
#
# TIPS ON USAGE
# 
# 1. This program should work on any system that can support perl.  The 
#    directory in which perl is located may vary by system.  On my system
#    it is /usr/bin/perl. The first line in this program should be altered 
#    to the correct directory for perl.
# 2. In UNIX, the program needs to be made executable by using 
#            chmod +x infocalc
# 3. In this program the -column flag indicates which column is used to 
#    identify individuals in the data file (either 1, 2, 3, 4, or 5).
#    The default is to use column 3.  Columns 1, 2, 4, and 5 are then 
#    ignored by the program.  However, it is necessary to have exactly 5
#    columns before the genotypes begin.
# 4. Negative numbers and 0 are treated as missing data, and are not
#    tabulated in calculating allele frequencies.
# 5. The program does not check for errors in formatting the data file
#    and may produce nonsensical results on incorrectly formatted data.  
#    However, since the format is the same as for STRUCTURE, you can run 
#    STRUCTURE in advance to check for errors.
# 6. The calculation of I_a assumes that there are at most 24 populations. 
# 7. The order in which the loci are printed in the output is alphabetical.
#
#########################################################################

###############################
###   End of instructions   ###
### The program begins here ###
###############################

###########################################
### Initialize and read in command-line ###
###########################################

$badconst = -9999;

use Getopt::Long;
my %opts = (column => 3, numpops => 2, input => "infile", output => "outfile", weightfile => "[none]");
Getopt::Long::GetOptions( \%opts, qw(
  column=i
  input=s
  numpops=i
  output=s
  weightfile=s
)) || die "Exiting... failed to parse command-line options\n";

if ($opts{column}<=0 || $opts{column}>5) 
  {die "Exiting... column specifying population identifiers (-column) must be in [1,5]\n";}
if ($opts{numpops}<=0) 
  {die "Exiting... number of populations (-numpops) must be a positive integer\n";}

$column  = $opts{column}-1;
$numPops = $opts{numpops};
$sample  = $opts{input};
$outfile = $opts{output};
$wfile   = $opts{weightfile};

open (OUTFI, ">$outfile") || die ("Exiting... can't open output file (-output) $outfile");

##########################################
### Read in file of population weights ###
##########################################

if ($wfile ne "[none]") {
  open (WEIGHTFI, $wfile) || die ("Exiting... can't open file of weights (-weightfile) $wfile");
  $sum=0;
  while ($line = <WEIGHTFI>) {
    chomp($line);
    @fields = split(/\s+/, $line);
    if ($#fields != 1) {die ("Exiting... weightfile $wfile does not have exactly two columns");}
    $pop = $fields[0];
    $weight{$pop} = $fields[1];
    $sum+=$weight{$pop};
  }
  if ( ($sum-1)*($sum-1) > 0.000001) {
    die ("Exiting... sum of weights in weightfile $wfile is not 1\n");
  }
}

#######################################################
### Read in list of loci from top line of data file ###
#######################################################

open (DATAFILE, $sample) || die ("Exiting... can't open input data file (-input) $sample");
$line = <DATAFILE>;
chomp($line);
@loci = split(/\s+/, $line);

##################################################################
### Read in data, obtain absolute frequencies and sample sizes ###
##################################################################

while ($line = <DATAFILE>) {
  chomp($line);
  @fields = split(/\s+/, $line);
  $pop = $fields[$column];
  if ($wfile eq "[none]") {$weight{$pop} = 1/$numPops;}
  $totalcount{$pop}+=0.5;
  for $i (5 .. $#fields) {
    $locus = $loci[$i-5];
    if ($fields[$i] > 0) {
      $allele = $fields[$i];
      ${${$absfreq{$locus}}{$allele}}{$pop}++;
      ${$sampsize{$locus}}{$pop}++;
    }
  }  
}

############################
### Perform a few checks ###
############################

$sum=0;
$num=0;
foreach $pop (keys %totalcount) {
  $sum+=$weight{$pop};
  if ($weight{$pop} < 0 || $weight{$pop} >= 1) {
    die ("Exiting... weight for population $pop is outside the range [0,1)\n");
  }
    $num++;
}  
if ($num != $numPops) {
  print ("Exiting... number of populations detected is not the same as in -numpops.\n");
  print ("  (1) Check that -numpops reflects the correct number of\n");
  print ("      populations in the appropriate column of the data file.\n");
  die   ("  (2) Check that the correct column was specified in -column.\n");
}
if ( ($sum-1)*($sum-1) > 0.000001) {
  print ("Exiting... the sum of population weights in the data file is not 1.\n");
  print ("  (1) Check that -numpops reflects the correct number of\n");
  print ("      populations in the appropriate column of the data file.\n");
  print ("  (2) Check that the correct column was specified in -column.\n");
  print ("  (3) If you are using a weightfile, check that the spelling \n");
  print ("      of population names is the same in the weightfile as in \n");
  die   ("      the data file.\n");
}
if ($numPops>24) { 
  {print "Warning: number of populations (-numpops) must be at most 24.\n";}
  {print "---> The computation of I_a will be omitted.\n";}
}

#####################################################################
### Compute relative frequencies for each locus-allele-population ###
### combination, and average frequencies across populations.      ###
#####################################################################

foreach $locus (keys %absfreq) {
  foreach $allele (keys %{$absfreq{$locus}}) {
    $aver = 0;
    foreach $pop (keys %{${$absfreq{$locus}}{$allele}}) {
      $r = ${${$absfreq{$locus}}{$allele}}{$pop}; 
      $s = ${$sampsize{$locus}}{$pop};
      ${${$relfreq{$locus}}{$allele}}{$pop} = $r / $s;
      $aver += $weight{$pop} * ${${$relfreq{$locus}}{$allele}}{$pop};
    }
    ${$averageRelfreq{$locus}}{$allele} = $aver;
  }
}

####################################################
### Compute informativeness for assignment (I_n) ###
####################################################

foreach $locus (keys %absfreq) {
  $left = 0;
  $right = 0;
  foreach $allele (keys %{$absfreq{$locus}}) {
    $pj = ${$averageRelfreq{$locus}}{$allele};
    if ($pj > 0) {$left -= $pj * log($pj)};
  }
  foreach $allele (keys %{$absfreq{$locus}}) {
    foreach $pop (keys %{${$absfreq{$locus}}{$allele}}) {   
      $pij = ${${$relfreq{$locus}}{$allele}}{$pop};
      if ($pij > 0) { $right += $weight{$pop} * $pij * log($pij); }
    }
  }
  $mutualinfo{$locus} = $left + $right;
}

############################################
### Compute ORCA (single-allele version) ###
############################################

foreach $locus (keys %absfreq) {
  $decisioninfo{$locus} = 0;
  foreach $allele (keys %{$absfreq{$locus}}) {
    $maxpijqi = 0;
    foreach $pop (keys %{${$absfreq{$locus}}{$allele}}) {   
      $pijqi = ${${$relfreq{$locus}}{$allele}}{$pop} * $weight{$pop};
      if ($pijqi > $maxpijqi) { $maxpijqi = $pijqi; }
    }
    $decisioninfo{$locus} += $maxpijqi;
  }
}

###############################################
### Compute ORCA (diploid genotype version) ###
###############################################

foreach $locus (keys %absfreq) {
  $dipldecisioninfo{$locus} = 0;
  foreach $allele (keys %{$absfreq{$locus}}) {
    foreach $otherallele (keys %{$absfreq{$locus}}) {
      $maxprod = 0;
      foreach $pop (keys %{${$absfreq{$locus}}{$allele}}) {   
        $prod  = ${${$relfreq{$locus}}{$allele}}{$pop} * $weight{$pop};
        $prod *= ${${$relfreq{$locus}}{$otherallele}}{$pop};
        if ($prod > $maxprod) { $maxprod = $prod; }
      }
      $dipldecisioninfo{$locus} += $maxprod;
    }
  }
}

##########################
### Compute factorials ###
##########################

BEGIN {
  my @fact = (1);
  sub factorial($) {
    my $n = shift;
    return $fact[$n] if defined $fact[$n];
    $fact[$n] = $n * factorial($n - 1);
  }
}

######################################
### List Stirling number constants ###
######################################

sub stirling {
  if ($numPops == 2) { return 3 };
  if ($numPops == 3) { return 11 };
  if ($numPops == 4) { return 50 };
  if ($numPops == 5) { return 274 };
  if ($numPops == 6) { return 1764 };
  if ($numPops == 7) { return 13068 };
  if ($numPops == 8) { return 109584 };
  if ($numPops == 9) { return 1026576 };
  if ($numPops == 10) { return 10628640 };
  if ($numPops == 11) { return 120543840 };
  if ($numPops == 12) { return 1486442880 };
  if ($numPops == 13) { return 19802759040 };
  if ($numPops == 14) { return 283465647360 };
  if ($numPops == 15) { return 4339163001600 };
  if ($numPops == 16) { return 70734282393600 };
  if ($numPops == 17) { return 1223405590579200 };
  if ($numPops == 18) { return 22376988058521600 };
  if ($numPops == 19) { return 431565146817638400 };
  if ($numPops == 20) { return 8752948036761600000 };
  if ($numPops == 21) { return 186244810780170240000 };
  if ($numPops == 22) { return 4148476779335454720000 };
  if ($numPops == 23) { return 96538966652493066240000 };
  if ($numPops == 24) { return 2342787216398718566400000 };
}

################################################
### Check if ok to compute I_a               ###
###   -need at most 24 populations           ###
###   -population weights need to be uniform ###
################################################

$computeIa=1;
if ($numPops > 24) {$computeIa=0;}
foreach $pop (keys %totalcount) {
  if ($weight{$pop} * $numPops != 1) {$computeIa=0;}
}

###############################################################
### Compute informativeness for ancestry coefficients (I_a) ###
###############################################################

if ($computeIa == 1) {
  $stirl = stirling($numPops);
  $facto = factorial($numPops);
  foreach $locus (keys %absfreq) {
    $sum = 0;
    foreach $allele (keys %{$absfreq{$locus}}) {
      $pj = ${$averageRelfreq{$locus}}{$allele};
      $sum += $pj * (1 - log($pj) - $stirl/$facto ); 
      foreach $pop (keys %{${$absfreq{$locus}}{$allele}}) {   
        $pij = ${${$relfreq{$locus}}{$allele}}{$pop};
        $num = 0;
        $denom = $numPops;
        if ($pij > 0) { 
          $num = ($pij ** $numPops) * log($pij);
          $z = 0;
          foreach $otherpop (keys %{${$absfreq{$locus}}{$allele}}) {   
            if ($otherpop ne $pop) {
              $plj = ${${$relfreq{$locus}}{$allele}}{$otherpop};
              $denom *= $pij - $plj;
            }
            $z++;
          }
          while ($z < $numPops) {
            $denom *= $pij;
            $z++;
          }
          if ($denom != 0 && $sum != $badconst) {
            $sum += $num/$denom;
          }
          if ($denom == 0) { 
            $sum = $badconst;
          }
        }
      }
    }
    $ancestryinfo{$locus} = $sum;
  }
}

############################
### Print the statistics ###
############################

printf OUTFI ("%10s\t%10s\t%10s\t", "Locus", "I_n", "I_a");
printf OUTFI ("ORCA[1-allele]\tORCA[2-allele]\n");
foreach $locus (sort keys %absfreq) {
  printf OUTFI ("%10s\t%10g\t", $locus, $mutualinfo{$locus});
  if ($computeIa==1) {printf OUTFI ("%10g\t", $ancestryinfo{$locus})};
  if ($computeIa==0) {printf OUTFI ("%10s\t", "NA")};
  printf OUTFI ("%10g\t", $decisioninfo{$locus});
  printf OUTFI ("%10g\n", $dipldecisioninfo{$locus});
}

########################################################
### Print summary of command line to the output file ###
########################################################

printf OUTFI "Command:       ";
printf OUTFI "infocalc  -column $opts{column}  ";
printf OUTFI "-numpops $opts{numpops}  ";
printf OUTFI "-input $opts{input}  ";
printf OUTFI "-output $opts{output}  ";
printf OUTFI "-weightfile $opts{weightfile}\n";

#################################################
### Print the weights used in the calculation ###
#################################################

printf OUTFI ("PriorWeights:  ");
foreach $pop (sort keys %totalcount) {
  printf OUTFI ("%s %g   ", $pop, $weight{$pop});
} printf OUTFI ("\n");


#################################################
#################################################
#################################################
#
# Notes
#
# [Version 1.1 incorporates new features from informativeStats-weights.perl.]
#
# [Version 1.0 was constructed from informativeStats2.perl for I_n and I_a,
# and from decisionTheory.perl for ORCA, removing nonessential features.] 

