#!/usr/bin/perl -w

use strict;

## We expect to be passed the 454 454PairStatus.txt file(s) for the
## Paired Ends.

die "pass 454PairStatus file\n"
  unless @ARGV;



## Read in the alignment results and build a list of dupes

my ($c, %dupe_list);

warn "\nparsing results\n";

while(<>){
  chomp;
  
  my ($clone, $status, $distance,
      $hn1, $hs1, $st1,
      $hn2, $hs2, $st2) = split/\t/;
  
  next unless
    $status eq "TruePair" ||
    $status eq "FalsePair";
  
  ## NB: Strandedness (st) never seems to be an issue with dupes,
  ## i.e. the same clone is amplified and resequenced, aligning with
  ## the same orientation. However, we keep strandedness in the check
  ## for completeness.
  
  my $uniq1 =
    "$hn1:$hs1:$st1".
    "$hn2:$hs2:$st2";
  
  $c++;
  push @{$dupe_list{$uniq1}}, $_;
}

warn "\ndiscovered ", scalar keys %dupe_list,
    " nr pairs out of $c total pairs\n";



warn "\ndumping\n";

for (keys %dupe_list){
  my $g = scalar @{$dupe_list{$_}};
  print @{$dupe_list{$_}}[0], "\t$g\n";
}

warn "OK\n";
