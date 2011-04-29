#!/usr/bin/perl -w

use strict;
use Getopt::Long;

=head1 DESCRIPTION

This script produces GFF from specific 454 PE alignment summaries. The
script doesn't do anything fancy. However, the analysis is based on MP
alignments in Newbler.

=cut



## OPTIONS

## Set the source and type fields, if desired
my $species = 'DM';
my $target  = 'PGSC0003DMB';
my $source  = 'dundee';
my $type    = 'match';

my $library = '';

## Do we have access to trim file(s)? (gives us trimmed read lengths)
my @trim_files;

## How about read file(s)? (gives us a read alignment summary)
my @read_files;

## Set to 1 to enable debugging output, if desired
my $verbose = 0;

GetOptions( "species=s"   => \$species,
            "version=s"   => \$target,
            "source=s"    => \$source,
            "type=s"      => \$type,
            
            "library=s"   => \$library,
            
            "trim=s{,}"   => \@trim_files,
            "read=s{,}"   => \@read_files,
            
            "verbose|v+"  => \$verbose,
          )
  or die "failed to communicate!\n";



## We expect to be passed the 454 454PairStatus.txt file(s) for the
## Paired Ends.

die "pass 454PairStatus file\n"
  unless @ARGV;



my %trim_data;
for(@trim_files){
  warn "Collecting trim data from $_\n"
    if $verbose > 0;
  open TF, '<', $_
    or die "fail : $!\n";
  
  while(<TF>){
    my @trim = split/\t/;
    ## Trimpoints Used / Used Trimmed Length
    $trim_data{$trim[0]} = [@trim[1,2]];
  }
}

warn "found ", scalar keys %trim_data, "\n"
  if $verbose > 0;



my %read_data;
for(@read_files){
  warn "Collecting read data from $_\n"
    if $verbose > 0;
  open RF, '<', $_
    or die "fail : $!\n";
  
  while(<RF>){
    my @read = split/\t/;
    ## Mapping Status / Mapped Accuracy / % of Read Mapped
    $read_data{$read[0]} = [@read[1,2,3]];
  }
}

warn "found ", scalar keys %read_data, "\n"
  if $verbose > 0;





## Read in the alignment results and produce the GFF

warn "\nparsing newbler results and dumping GFF\n";

while(<>){
  
  my ($clone, $status, $distance,
      $hit_name1, $hit_beg1, $strand1,
      $hit_name2, $hit_beg2, $strand2) = split/\s+/;
  
  # debugging
  die "failed to find trim or read data for $clone\n"
    unless
      exists $trim_data{"$clone\_left" } &&
      exists $trim_data{"$clone\_right"} &&
      exists $read_data{"$clone\_left" } &&
      exists $read_data{"$clone\_right"};
  
  
  
  ## Get additional read mapping data
  my ($map_status1, $map_accuracy1, $perc_mapped1) = @{$read_data{"$clone\_left" }};
  my ($map_status2, $map_accuracy2, $perc_mapped2) = @{$read_data{"$clone\_right"}};
  
  ## Only want 'Full' matches!
  next if $map_status1 ne "Full";
  next if $map_status2 ne "Full";
  
  
  ## Get read trim data
  my ($trim_points1, $trim_length1) = @{$trim_data{"$clone\_left" }};
  my ($trim_points2, $trim_length2) = @{$trim_data{"$clone\_right" }};
  
  ## Read end points are not recorded in the pair file!
  my $hit_end1 = $hit_beg1 + $trim_length1;
  my $hit_end2 = $hit_beg2 + $trim_length2;
  
  my ($query_beg1, $query_end1) = split(/-/, $trim_points1);
  my ($query_beg2, $query_end2) = split(/-/, $trim_points2);
  
  my $score1 = $trim_length1 * $map_accuracy1 / 100;
  my $score2 = $trim_length2 * $map_accuracy2 / 100;
  
  print
    join("\t",
         $hit_name1, $source, $type, $hit_beg1, $hit_end1, $score1, $strand1, '.',
         join(";", "ID=$library$clone\_F", "Target=x $query_beg1 $query_end1", "ident=$map_accuracy1")
        ), "\n";
  
  print
    join("\t",
         $hit_name2, $source, $type, $hit_beg2, $hit_end2, $score2, $strand2, '.',
         join(";", "ID=$library$clone\_R", "Target=x $query_beg2 $query_end2", "ident=$map_accuracy2")
        ), "\n";
  
}

warn "OK\n";
