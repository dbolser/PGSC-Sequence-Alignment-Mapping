#!/usr/bin/perl -w

use strict;

die "foo you!\n"
  unless @ARGV == 4;

my $SSAHA_R  = $ARGV[0];
my $SSAHAIDX = $ARGV[1];
my $QUERY    = $ARGV[2];
my $RESULTS  = $ARGV[3];

open O, ">", "$RESULTS"
  or die "fooie!\n";

open P, "-|", "$SSAHA_R -tags 0 -save $SSAHAIDX $QUERY"
  or die "fooie2!\n";

die "$!\n" if $?;

while(<P>){
  if (/^\d/){
    s/ +/\t/g;
    print O;
  }
  elsif (/^SSAHA2:/){
    print;
  }
  else{
    ## Other lines are not important
  }
}

warn "OK\n";
