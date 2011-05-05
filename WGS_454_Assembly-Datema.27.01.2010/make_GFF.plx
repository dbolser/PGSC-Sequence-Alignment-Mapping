#!/usr/bin/perl -w

## Script to convert an ATAC A2A MAPPER clump file into GFF format

use strict;
use Getopt::Long;

my $source;

my $verbose = 0;
my $type    = 'match';

GetOptions(
           "verbose"    => \$verbose,
           "source|s=s" => \$source,
           "type|t=s"   => \$type,
          )
  or die "failed to parse command line options\n";

die "set the GFF source field with --source '<source>'\n"
  unless defined($source);

my $clump_file = $ARGV[0];

die "pass a clump file you lump\n"
  unless -s $clump_file;





## Go time

my $rows;

while(<>){
  unless(/^M/){
    warn "skipping header line : $_"
      if $verbose > 0;
    next
  }
  
  ## Grab the clumps
  next unless /^M c /;
  
  $rows++;
  
  my ($class, $match_type, $match_id, $pmatch_id,
      $hit, $h_off, $h_len, $h_strand,
      $que, $q_off, $q_len, $q_strand,
     ) = split /\s/;
  
  my $h_end = $h_off + $h_len;
  my $q_end = $q_off + $q_len;
  
  ## Pick 'direction' of GFF...
  
#   print
#     join("\t",
#            $que, $source, $type, $q_off, $q_end, 0, '.', ($q_strand > 0 ? '+' : '-'),
#          join(";",
#               "Name=$hit",
#               "Target=$hit $h_off $h_end $h_strand"
#              ),
#         ), "\n";
  
  print
    join("\t",
         $hit, $source, $type, $h_off, $h_end, 0, '.', ($q_strand > 0 ? '+' : '-'),
         join(";",
              "ID=$que",
              "Name=$que",
              "Target=$que $q_off $q_end $q_strand"
             ),
        ), "\n";

}

warn "$rows rows of GFF\n";
warn "OK\n";

