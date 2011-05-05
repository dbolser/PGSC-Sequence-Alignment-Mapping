#!/usr/bin/perl -w

use strict;

my $verbose = 0;

my $clump_file = $ARGV[0];

die "pass a clump file you lump\n"
  unless -s $clump_file;

## output
open O1, '>', 'out1.csv'
  or die "failed to open 'out1.csv' for writing\n";
open O2, '>', 'out2.csv'
  or die "failed to open 'out2.csv' for writing\n";



## Grab sequence lengths

my $psl_file = 'Data/Potato_Assembly/PGSC0003DM/Indexes/PGSC0003DMB.fa.len';
my $tsl_file = 'Assembly/S_lycopersicum_chromosomes.2.10.fa.len';


open PSL, '<', $psl_file or die "failed to open $psl_file : $!\n";
open TSL, '<', $tsl_file or die "failed to open $tsl_file : $!\n";

my %psl;
while(<PSL>){chomp; my ($id, $len) = split/\t/; $psl{$id} = $len}
warn "got the length of ", scalar keys %psl, " potato sequences\n";

my %tsl;
while(<TSL>){chomp; my ($id, $len) = split/\t/; $tsl{$id} = $len}
warn "got the length of ", scalar keys %tsl, " tomato sequences\n";



## Parse clumpy

my @hits;

while(<>){
  unless(/^M/){
    warn "skipping header line : $_"
      if $verbose > 0;
    next
  }
  
  ## Grab the clumps
  next unless /^M c /;
  
  my ($class, $match_type, $match_id, $pmatch_id,
      $hit, $h_off, $h_len, $h_strand,
      $que, $q_off, $q_len, $q_strand,
     ) = split /\s/;
  
  die "seq?\t$hit\n" unless my $hs_len = $tsl{$hit};
  die "seq?\t$que\n" unless my $qs_len = $psl{$que};
  
  
  my $h_end = $h_off + $h_len;
  my $q_end = $q_off + $q_len;
  
  my $h_offx = sprintf("%06.3f", $h_off/$hs_len * 100);
  my $h_endx = sprintf("%06.3f", $h_end/$hs_len * 100);
  my $h_lenx = sprintf("%06.3f", $h_len/$hs_len * 100);
  
  my $q_offx = sprintf("%06.3f", $q_off/$qs_len * 100);
  my $q_endx = sprintf("%06.3f", $q_end/$qs_len * 100);
  my $q_lenx = sprintf("%06.3f", $q_len/$qs_len * 100);
  
  push @hits,
    [$hit, $hs_len, $h_off, $h_end, $h_len, $h_offx, $h_endx, $h_endx - $h_offx, $h_strand, '',
     $que, $qs_len, $q_off, $q_end, $q_len, $q_offx, $q_endx, $q_endx - $q_offx, $q_strand
    ];
  
}

warn "OK\n";



my $tog1 = '';
for my $hit (sort sit @hits){
  print O1 "\n" if $hit->[0] ne $tog1; $tog1 = $hit->[0];
  print O1
    join("\t", @{$hit}[0..8], '', @{$hit}[10..18]), "\n";
}
warn "OK\n";

my $tog2 = '';
for my $hit (sort sat @hits){
  print O2 "\n" if $hit->[10] ne $tog2; $tog2 = $hit->[10];
  print O2
    join("\t", @{$hit}[10..18], '', @{$hit}[0..8]), "\n";
}
warn "OK\n";



sub sit {
  return
    $a->[0] cmp $b->[0] ||
    $a->[2] <=> $b->[2] ||
    $b->[3] <=> $a->[3] ||
      (warn "eep\n" && 0);
}



sub sat {
  return
    $a->[10] cmp $b->[10] ||
    $a->[12] <=> $b->[12] ||
    $b->[13] <=> $a->[13] ||
      (warn "eep\n" && 0);
}
