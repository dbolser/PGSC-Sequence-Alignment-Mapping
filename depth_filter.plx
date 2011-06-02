#!/usr/bin/perl -w

## A script to check the 'depth' of alignment at every position on
## every reference(hit) and query sequence in a blast tabular
## report. Alignments within regions above the 'depth threshold' are
## discarded.

## NB: blastable format or ssaha format is hard-coded below.



use strict;
use Getopt::Long;

## DEFAULT: Regions of alignments with a 'depth' this deep (or deeper)
## are deemed 'repetative'.
my $depthThreshold = 3;

## DEFAULT: The fraction of the alignment that must have a depth BELOW
## the depth threshold to be acceptable.
my $coverageThreshold = 0.8;

my $file;

my $verbose = 0;

## Parse the command line for options
GetOptions(
           "depth=i"    => \$depthThreshold,
           "coverage=f" => \$coverageThreshold,
           "infile=s"   => \$file,
           "verbose"    => \$verbose,
          )
  or die "could not parse command line for options\n";

$file = $ARGV[0]
  unless $file;

die "depth should be a positive integer greater than 1\n"
  unless $depthThreshold > 1;

die "coverage should be a fraction in the range 0.10 to 0.90\n"
  unless $coverageThreshold > 0 && $coverageThreshold < 1;

die "give me an infile!\n"
  unless -s $file;



warn "using the following settings:\n";
warn "depth threshold : $depthThreshold\n";
warn "coverage threshold : $coverageThreshold\n";
warn "\n";



## Go time

## Calculate the positional 'depth' accross this sequence

my %depth;

warn "parsing alignment report\n";

open DATA, '<', $file
  or die "failed to open input file : $!\n";

while(<DATA>){
  my ($q, undef, undef, $h1, $p1, undef, $h2, $p2, undef)
    = split/\t/;
  
  ## I can't work it out from the manual... A picture would be nice!
  my $hs1 = $p1 - 200;
  my $he1 = $p1 + 200;
  
  my $hs2 = $p2 - 200;
  my $he2 = $p2 + 200;
  
  $hs1 = $hs1 < 0 ? 0 : $hs1;
  $hs2 = $hs2 < 0 ? 0 : $hs2;
  
  for my $i ( $hs1 .. $he1 ){
    $depth{$h1}[$i]++
  }
  for my $i ( $hs2 .. $he2 ){
    $depth{$h2}[$i]++
  }
}

warn "created depth vector for ", scalar( keys %depth ), " sequencs\n";



warn "applying the filter\n";

open DATA, '<', $file
  or die "failed to open input file : $!\n";

while(<DATA>){
  my ($q, undef, undef, $h1, $p1, undef, $h2, $p2, undef)
    = split/\t/;
  
  ## I can't work it out from the manual... A picture would be nice!
  my $hs1 = $p1 - 200;
  my $he1 = $p1 + 200;
  
  my $hs2 = $p2 - 200;
  my $he2 = $p2 + 200;
  
  $hs1 = $hs1 < 0 ? 0 : $hs1;
  $hs2 = $hs2 < 0 ? 0 : $hs2;
  
  # Filter alignments by depth!
  
  my ($hl, $ht);
  for my $i ($hs1 .. $he1){
    $hl++;
    $ht++ if $depth{$h1}[$i] >= $depthThreshold;
  }
  for my $i ($hs2 .. $he2){
    $hl++;
    $ht++ if $depth{$h2}[$i] >= $depthThreshold;
  }
  
  #next if defined($ht); # Simulate $coverateThrshold = 0;
  next if ($ht || 0) / $hl > ( 1 - $coverageThreshold );
  
  # We passed the filter...
  print "$_";
}

warn "OK\n";
