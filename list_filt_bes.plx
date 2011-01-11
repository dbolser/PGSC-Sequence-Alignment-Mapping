#!/usr/bin/perl -w

use strict;


my $dir = "../Results";

my $file = "PGSC0003DMB_vs_Solanum_phureja_DM.bac_ends.09242009.ssaha";

my $filt = "dfilt.q.d5.c80";

open L, '<', "$dir/$file"
  or die "fail 1\n";

open T, '<', "$dir/$file.$filt"
  or die "fail 2\n";

warn "parsing full\n";
my %full;
while(<L>){
  my $bes = (split)[1];
  die "'$bes'\n" unless $bes =~ /^gi\|(\d+)\|gb\|GS\d+\.1\|$/;
  $full{$1}++;
}
warn "got ", scalar keys %full, " BES full\n";


warn "parsing filt\n";
my %filt;
while(<T>){
  my $bes = (split)[1];
  die "'$bes'\n" unless $bes =~ /^gi\|(\d+)\|gb\|GS\d+\.1\|$/;
 $filt{$1}++;
}
warn "got ", scalar keys %filt, " BES filt\n";


warn "dumping the missing BES\n";

for my $bes_id (keys %full){
  print "$bes_id\n"
    unless exists $filt{$bes_id};
}
