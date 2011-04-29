#!/usr/bin/perl -w

use strict;
use Getopt::Long;

=head1 DESCRIPTION

This script produces GFF from SSAHA2 format alignments. The script
doesn't do anything fancy. Could I just use the standard BioPerl
script to do this?

=cut



## OPTIONS

## Set the source and type fields, if desired
my $species = 'DM';
my $target  = 'PGSC0003DMB';
my $source  = 'dundee';
my $type    = 'match';

## Define a 'query name mapping' file, usually derived from the fasta
## headerline for the query sequences, if desired
my @header_file;

## Set to 1 to enable debugging output, if desired
my $verbose = 0;

GetOptions( "species=s"   => \$species,
            "version=s"   => \$target,
            "source=s"    => \$source,
            "type=s"      => \$type,
            
            "header=s{3}" => \@header_file,
            
            "verbose|v+"  => \$verbose,
          )
  or die "failed to communicate!\n";



## We expect to be passed the (filtered) SSAHA2 alignment file(s)

die "pass a ssaha file\n"
  unless @ARGV;





## header file processing

my ($header_file,
    $header_file_match_col,
    $header_file_alias_col
   ) = @header_file;

my %mapping;

if(@header_file){
  die "cant open header file '$header_file' : $!\n"
    unless -s $header_file;
  
  unless( $header_file_match_col eq int $header_file_match_col  and
          $header_file_alias_col eq int $header_file_alias_col ){
    die "need integers for these values : ",
      "'$header_file_match_col' and '$header_file_alias_col'!\n";
  }
  
  open M, '<', $header_file
    or die "failed to open file '$header_file' : $!\n";
  
  while(<M>){
    chomp;
    next if /^#/ || /^\s*$/;
    
    my @row = split/\t/;
    
    die "redundant mapping!\n"
      if exists $mapping{$row[$header_file_match_col]};
    
    $mapping{$row[$header_file_match_col]}
      = $row[$header_file_alias_col];
  }
  
  warn "got ", scalar keys %mapping,
    " mapped headers from '$header_file'\n";
}





## Read in the alignment results and produce the GFF

warn "\nparsing ssaha results and dumping GFF\n";

print "##gff-version 3\n";
print "##species $species\n";
print "##genome-build $target\n";



my $hsp_counter;

while(<>){
  my ($score, $query_name, $hit_name, $qs, $qe, $hs, $he,
      $strand, $align_len, $ident, $query_len) = split/\s+/;
  
  $hsp_counter++;
  
  my $attrs;
  
  if(@header_file){
    my $query_name2 = $mapping{$query_name}
      or die "FAILED! : '$query_name'\n";
    
    $attrs =
      join(';',
           "ID=$query_name2",
           "Alias=$query_name",
           "Target=$query_name2 $qs $qe",
          );
  }
  else{
    $attrs =
      join(';',
           "ID=$query_name",
           "Target=$query_name $qs $qe",
          );
  }
  
  $strand =~ tr/FC/+-/
    or die "unexpected 'strand' : '$strand'\n";
  
  print
    join("\t",
         $hit_name,
         $source,
         $type,
         $hs, $he,
         $score,
         $strand,
         '.',
         $attrs. "ident=$ident;length=$query_len",
        ), "\n";
}

warn "got $hsp_counter HSPs\n";
warn "OK\n";
