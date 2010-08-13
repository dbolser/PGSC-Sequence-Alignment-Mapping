#!/usr/bin/perl -w

use strict;

=head1 DESCRIPTION

Thsi script produces GFF from specific BES or FES alignments. It uses
a specific, hard-coded regexp to match one end to another.

=cut



## Get an important file

use Getopt::Long;

## Set to 1 to enable debugging output.
my $verbose = 0;

## Table used to link forward and reverse sequence IDs to their common
## clone ID i.e. for paring forward and reverse sequences.
my $header_table;

GetOptions(
	   "v"    => \$verbose,
	   "hf=s" => \$header_table,
	  )
  or die "failed to parse command line options!\n";


die "pass a header table\n" unless defined($header_table);


## We expect to be passed the SSAHA2 alignment file for the FES
die "pass ssaha file\n"
  unless @ARGV;



## Read in the sequence header table so that we can match the sequence
## identifier in the alignment file to the more informative sequence
## header (parsed into the table elsewhere).

## At this stage, we also construct the FES pairs for a Fosmid (clone)

warn "\nreading sequence data table ($header_table)\n";

open SEQ, '<', $header_table
  or die "failed to open $header_table : $!\n";

my (%sequence_pairs, $seq_count);

while(<SEQ>){
  my ($gi, $gb, $gb_ver, $clone,
      $direction, $length) = split/\t/;
  
  $seq_count++;
  
  ## There should only be two directions!
  die unless
    $direction eq 'TF' ||
    $direction eq 'TR';
  
  ## No direction should be duplicated!
  die if
    exists($sequence_pairs{$clone}{$direction});
  
  ## Store the 'putative pair'
  $sequence_pairs{$clone}{$direction} = $gi;
}

warn "got $seq_count FES (and ",
  scalar keys %sequence_pairs, " Fosmids)\n";

warn "\nNote, not all Fosmids have FES for both directions!\n";





## Read in the alignment results

warn "\nparsing ssaha results\n";

my (%hits, $hit_counter);

while(<>){
  my ($score, $query_name, $hit_name, $qs, $qe, $hs, $he,
      $strand, $align_len, $ident, $query_len) = split/\s+/;
  
  $hit_counter++;
  
  ## Get the GI
  die "failed to parse : '$query_name'\n"
    unless $query_name =~ /^gi\|(\d+)\|gb\|(FI\d+)\.\d\|\2$/;
  
  ## Store all the hits for this GI together
  push
    @{$hits{$1}},
      [$score, $query_name, $hit_name, $qs, $qe, $hs, $he,
       $strand, $align_len, $ident, $query_len, $2];
}

warn "got $hit_counter hits for ",
  scalar keys %hits, " FES\n";





## Produce the GFF, paired and all

warn "\ndumping GFF\n";


CLONE:
foreach my $clone (keys %sequence_pairs){
  
  ## Really debugging
  #next unless $clone eq 'LuSp113E12';
  
  ## Debugging
  print "C:\t$clone\n"
    if $verbose > 0;
  
  ## Get the hits for both sequences
  my $gi1 = $sequence_pairs{$clone}{'TF'} || 'seq missing';
  my $gi2 = $sequence_pairs{$clone}{'TR'} || 'seq missing';
  
  my @hits1 = exists($hits{$gi1}) ? @{$hits{$gi1}} : ();
  my @hits2 = exists($hits{$gi2}) ? @{$hits{$gi2}} : ();
  
  ## Debugging
  print "G:\t$gi1\t", scalar(@hits1), "\n"
    if $verbose > 0;
  print "G:\t$gi2\t", scalar(@hits2), "\n"
    if $verbose > 0;
  
  
  
  ## Create some failure flags?
  my ($missing_f, $ali_f) = (0, 0);
  
  if ($gi1 eq 'seq missing' ||
      $gi2 eq 'seq missing'){
    $missing_f++;
  }
  
  if (@hits1 == 0){
    $ali_f++;
  }
  if (@hits2 == 0){
    $ali_f++;
  }
  
  ## OK
  
  
  
  
  
  ## Try to pair the hits (best hits first)
  
  ## NOTE: if the best hits don't form a 'good pair', we give up on
  ## the clone.

  my ($scaff_f, $ori_f, $short_f, $long_f) = (0, 0, 0, 0);
  
 HIT:
  for my $hit1 (sort hits_by_score @hits1){
    for my $hit2 (sort hits_by_score @hits2){
      
      my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1,
	  $al1, $nt1, $ql1, $gb1) = @$hit1;
      my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2,
	  $al2, $nt2, $ql2, $gb2) = @$hit2;
      
      ## Debugging
      print join("\t", '1:', @$hit1), "\n"
	if $verbose > 0;
      print join("\t", '2:', @$hit2), "\n"
	if $verbose > 0;
      
      ## Check they hit the same scaffold
      if ($hn1 ne $hn2){
	## Debugging
	print "W:\tnot on the same scaffold\n"
	  if $verbose > 0;
	$scaff_f++;
	last HIT;
      }
      
      ## Check orientation 1/2
      if ($st1 eq 'F'){
	if ($st2 ne 'C'){
	  ## Debugging
	  print "W:\tincorrect orientation\n"
	    if $verbose > 0;
	  $ori_f++;
	  last HIT;
	}
	
	## OK. Check their distance
	my $fes_distance = $he2 - $hs1 + 1;
	print "I:\tdistance is $fes_distance\n"
	  if $verbose > 0;
	
	if ($fes_distance < 5000){
	  ## Debugging
	  print "W:\treally incorrect orientation\n"
	    if $verbose > 0;
	  $short_f++;
	  last HIT;
	}
	
	if ($fes_distance > 100_000){
	  ## Debugging
	  print "W:\ttoo long!\n"
	    if $verbose > 0;
	  $long_f++;
	  last HIT;
	}
	
	## OK! Print three lines of GFF
	print join("\t", $hn1, 'dundee', 'Fosmid', $hs1, $he2,  '.', '+', '.', "ID=$clone;Name=$clone"), "\n";
	print join("\t", $hn1, 'dundee', 'FES',    $hs1, $he1, $sc1, '+', '.', "ID=$gi1;Name=$gb1;Parent=$clone;Note=TF"), "\n";
	print join("\t", $hn1, 'dundee', 'FES',    $hs2, $he2, $sc2, '-', '.', "ID=$gi2;Name=$gb2;Parent=$clone;Note=TR"), "\n";
	
	## DONE!
	next CLONE;
      }
      
      
      
      ## Check orientation 2/2
      if ($st1 eq 'C'){
	if ($st2 ne 'F'){
	  ## Debugging
	  print "W:\tincorrect orientation\n"
	    if $verbose > 0;
	  $ori_f++;
	  last HIT;
	}
	
	## OK. Check their distance
	my $fes_distance = $he1 - $hs2 + 1;
	print "I:\tdistance is $fes_distance\n"
	  if $verbose > 0;
	
	if ($fes_distance < 5000){
	  ## Debugging
	  print "W:\treally incorrect orientation\n"
	    if $verbose > 0;
	  $short_f++;
	  last HIT;
	}
	
	if ($fes_distance > 100_000){
	  ## Debugging
	  print "W:\ttoo long!\n"
	    if $verbose > 0;
	  $long_f++;
	  last HIT;
	}
	
	## OK! Print three lines of GFF
	print join("\t", $hn1, 'dundee', 'Fosmid', $hs2, $he1,  '.', '-', '.', "ID=$clone;Name=$clone"), "\n";
	print join("\t", $hn1, 'dundee', 'FES',    $hs1, $he1, $sc1, '-', '.', "ID=$gi1;Name=$gb1;Parent=$clone;Note=TF"), "\n";
	print join("\t", $hn1, 'dundee', 'FES',    $hs2, $he2, $sc2, '+', '.', "ID=$gi2;Name=$gb2;Parent=$clone;Note=TR"), "\n";
	
	## DONE!
	next CLONE;
      }
      
    }
  }
  
  
  
  
  
  
  
  ## Nothing paired, falling through to single end printer\n";
  
  ## NB! We simply don't care about a whole class of failed pairs!
  next if $missing_f || $ali_f;
  
  ## NB! We simply don't care about another class of failed pairs!
  next if $ori_f || $short_f || $long_f;
  
  
  
  ## The only remaining unpaired ends are those that span
  ## superscaffolds!
  
  ## Remeber 1 = TF = forward
  ## Remeber 2 = TR = reverse
  my $hit1 = (sort hits_by_score @hits1)[0];
  my $hit2 = (sort hits_by_score @hits2)[0];
  
  my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1, $al1, $nt1, $ql1, $gb1) = @$hit1;
  my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2, $al2, $nt2, $ql2, $gb2) = @$hit2;
  
  ## Debugging
  print join("\t", '1:', @$hit1), "\n" if $verbose > 0;
  print join("\t", '2:', @$hit2), "\n" if $verbose > 0;
  
  ## ONE
  if($st1 eq 'F'){
    print join("\t", $hn1, 'dundeex', 'Fosmid', $hs1, $he1+1000,  '.', '+', '.',
	       "ID=$clone.TF;Name=$clone.TF;Note=Other end matches $hn2"), "\n";
    print join("\t", $hn1, 'dundeex', 'FES',    $hs1, $he1,      $sc1, '+', '.',
	       "ID=$gi1;Name=$gb1;Parent=$clone.TF"), "\n";
  }
  if($st1 eq 'C'){
    print join("\t", $hn1, 'dundeex', 'Fosmid', $hs1-1000, $he1,  '.', '-', '.',
	       "ID=$clone.TF;Name=$clone.TF;Note=Other end matches $hn2"), "\n";
    print join("\t", $hn1, 'dundeex', 'FES',    $hs1, $he1,      $sc1, '-', '.',
	       "ID=$gi1;Name=$gb1;Parent=$clone.TF"), "\n";
  }
  
  ## TWO
  if($st2 eq 'F'){
    print join("\t", $hn2, 'dundeex', 'Fosmid', $hs2, $he2+1000,  '.', '-', '.',
	       "ID=$clone.TR;Name=$clone.TR;Note=Other end matches $hn1"), "\n";
    print join("\t", $hn2, 'dundeex', 'FES',    $hs2, $he2,      $sc2, '+', '.',
	       "ID=$gi2;Name=$gb2;Parent=$clone.TR"), "\n";
  }
  if($st2 eq 'C'){
    print join("\t", $hn2, 'dundeex', 'Fosmid', $hs2-1000, $he2,  '.', '+', '.',
	       "ID=$clone.TR;Name=$clone.TR;Note=Other end matches $hn1"), "\n";
    print join("\t", $hn2, 'dundeex', 'FES',    $hs2, $he2,      $sc2, '-', '.',
	       "ID=$gi2;Name=$gb2;Parent=$clone.TR"), "\n";
  }
  
}

warn "OK\n";



sub hits_by_score {
  # score
  $b->[0] <=> $a->[0] ||
    # ident
    $b->[9] <=> $a->[9];
}

