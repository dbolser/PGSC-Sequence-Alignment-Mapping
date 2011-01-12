#!/usr/bin/perl -w

use strict;

=head1 DESCRIPTION

This script produces GFF from BES alignments (SSAHA2 format). It uses
a specific regexp to match one end to another, and uses a BES sequence
specific data table to construct the GFF. The script removes very
short or very long pairs, so they will never be seen!

=cut


## We expect to be passed the SSAHA2 alignment file for the BES
die "pass ssaha file\n"
  unless @ARGV;



## Set to 1 to enable debugging output
my $verbose = 0;



## Read in the BES specif sequence data table so that we can match the
## sequence identifier in the alignment file to the more informative
## sequence header (the header is parsed into a table elsewhere).

## At this stage, we also construct the BES pairs for a BAC (clone)

warn "\nreading sequence data table\n";

my $seq_data_file =
  "Data/Solanum_phureja_DM.bac_ends.09242009.header.tab";

open SEQ, '<', $seq_data_file
  or die "failed to open $seq_data_file : $!\n";

my (%sequence_pairs, $seq_count);

while(<SEQ>){
  my ($gi, $gb, $gb_ver, $clone, $clone_id,
      $direction, $length) = split/\t/;
  
  $seq_count++;
  
  ## There should only be two directions!
  die unless
    $direction eq 'TP' ||
    $direction eq 'TV';
  
  ## No direction should be duplicated!
  die if
    exists($sequence_pairs{$clone}{$direction});
  
  ## Store the 'putative pair'
  $sequence_pairs{$clone}{$direction} = $gi;
}

warn "got $seq_count BES (and ",
  scalar keys %sequence_pairs, " BACs)\n";

warn "\nNote, not BACs have BES for both directions!\n";





## Read in the alignment results

warn "\nparsing ssaha results\n";

my (%hits, $hit_counter);

while(<>){
  my ($score, $query_name, $hit_name, $qs, $qe, $hs, $he,
      $strand, $align_len, $ident, $query_len) = split/\s+/;
  
  $hit_counter++;
  
  ## Get the GI
  die unless $query_name =~ /^gi\|(\d+)\|gb\|GS\d+\.1\|$/;
  
  ## Store all the hits for this GI together
  push
    @{$hits{$1}},
      [$score, $query_name, $hit_name, $qs, $qe, $hs, $he,
       $strand, $align_len, $ident, $query_len];
}

warn "got $hit_counter hits for ",
  scalar keys %hits, " BES\n";







## Produce the GFF

warn "\ndumping GFF\n";

foreach my $clone (keys %sequence_pairs){
  
  ## Really debugging
  #next unless $clone eq 'LuSp113E12';
  
  ## Debugging
  print "C:\t$clone\n"
    if $verbose > 0;
  
  ## Get the hits for both sequences
  my $gi1 = $sequence_pairs{$clone}{'TP'} || 'seq missing';
  my $gi2 = $sequence_pairs{$clone}{'TV'} || 'seq missing';
  
  my @hits1 = exists($hits{$gi1}) ? @{$hits{$gi1}} : ();
  my @hits2 = exists($hits{$gi2}) ? @{$hits{$gi2}} : ();
  
  ## Debugging
  print "G:\t$gi1\t", scalar(@hits1), "\n"
    if $verbose > 0;
  print "G:\t$gi2\t", scalar(@hits2), "\n"
    if $verbose > 0;
  
  ## Hits for both ends?
  next unless exists($hits{$gi1});
  next unless exists($hits{$gi2});
  
  
  
  
  
  ## Try to pair the hits (best hits first)
  
  ## NOTE: if the best hits don't form a 'good pair', we give up on
  ## the clone.
  
  ## NOTE: about 3% of the clones could be 'recovered' by hunting for
  ## a lower scoring 'good pair', however, their quality is dubious,
  ## and they are complex to handle! We skip them here.
  
  my $hit1 = (sort hits_by_score @hits1)[0]; # Remeber 1 = TP = forward
  my $hit2 = (sort hits_by_score @hits2)[0]; # Remeber 2 = TV = reverse
  
  my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1,
      $al1, $nt1, $ql1) = @$hit1;
  my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2,
      $al2, $nt2, $ql2) = @$hit2;
  
  ## Debugging
  print join("\t", '1:', @$hit1), "\n"
    if $verbose > 0;
  print join("\t", '2:', @$hit2), "\n"
    if $verbose > 0;
  
  
  
  ## Check they hit the same scaffold (if not they must span a pair of
  ## scaffolds)
  
  if ($hn1 eq $hn2){
    
    ## Check orientation 1/2
    if ($st1 eq 'F'){
      if ($st2 ne 'C'){
	## Debugging
	print "W:\tincorrect orientation\n"
	  if $verbose > 0;
	next;
      }
      
      ## OK. Check their distance
      my $bes_distance = $he2 - $hs1 + 1;
      print "I:\tdistance is $bes_distance\n"
	if $verbose > 0;
      
      if ($bes_distance < 5_000){
	## Debugging
	print "W:\tincorrect orientation\n"
	  if $verbose > 0;
	  next;
      }
      
      if ($bes_distance > 500_000){
	## Debugging
	print "W:\ttoo long!\n"
	  if $verbose > 0;
	next;
      }
      
      ## OK! Print three lines of GFF
      print join("\t", $hn1, 'dundee', 'BAC', $hs1, $he2,  '.', '+', '.', "ID=$clone;Name=$clone"), "\n";
      print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1, $sc1, '+', '.', "ID=$gi1;Parent=$clone;Note=TP"), "\n";
      print join("\t", $hn1, 'dundee', 'BES', $hs2, $he2, $sc2, '-', '.', "ID=$gi2;Parent=$clone;Note=TV"), "\n";
      
      ## DONE!
      next;
    }
    
    
    
    ## Check orientation 2/2
    if ($st1 eq 'C'){
      if ($st2 ne 'F'){
	## Debugging
	print "W:\tincorrect orientation\n"
	  if $verbose > 0;
	next;
      }
      
      ## OK. Check their distance
      my $bes_distance = $he1 - $hs2 + 1;
      print "I:\tdistance is $bes_distance\n"
	if $verbose > 0;
      
      if ($bes_distance < 5_000){
	## Debugging
	print "W:\tincorrect orientation\n"
	  if $verbose > 0;
	next;
      }
      
      if ($bes_distance > 500_000){
	## Debugging
	print "W:\ttoo long!\n"
	  if $verbose > 0;
	next;
      }
      
      ## OK! Print three lines of GFF
      print join("\t", $hn1, 'dundee', 'BAC', $hs2, $he1,  '.', '-', '.', "ID=$clone;Name=$clone"), "\n";
      print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1, $sc1, '-', '.', "ID=$gi1;Parent=$clone;Note=TP"), "\n";
      print join("\t", $hn1, 'dundee', 'BES', $hs2, $he2, $sc2, '+', '.', "ID=$gi2;Parent=$clone;Note=TV"), "\n";
      
      ## DONE!
      next;
    }
  }
  
  
  
  ## Nothing paired, falling through to single end printer
  
  ## Debugging
  print "W:\tnot on the same scaffold\n"
    if $verbose > 0;
  
  ## The only remaining unpaired ends are those that span
  ## superscaffolds!
  
  ## Debugging
  print join("\t", '1:', @$hit1), "\n" if $verbose > 0;
  print join("\t", '2:', @$hit2), "\n" if $verbose > 0;
  
  ## ONE
  if($st1 eq 'F'){
    print join("\t", $hn1, 'dundeex', 'BAC', $hs1, $he1+1000,  '.', '+', '.',
	       "ID=$clone.TP;Name=$clone.TP;Note=Other end matches $hn2"), "\n";
    print join("\t", $hn1, 'dundeex', 'BES', $hs1, $he1,      $sc1, '+', '.',
	       "ID=$gi1;Parent=$clone.TP;Note=TP"), "\n";
  }
  if($st1 eq 'C'){
    print join("\t", $hn1, 'dundeex', 'BAC', $hs1-1000, $he1,  '.', '-', '.',
	       "ID=$clone.TP;Name=$clone.TP;Note=Other end matches $hn2"), "\n";
    print join("\t", $hn1, 'dundeex', 'BES', $hs1, $he1,      $sc1, '-', '.', 
	       "ID=$gi1;Parent=$clone.TP;Note=TP"), "\n";
  }
  
  ## TWO
  if($st2 eq 'F'){
    print join("\t", $hn2, 'dundeex', 'BAC', $hs2, $he2+1000,  '.', '-', '.',
	       "ID=$clone.TV;Name=$clone.TV;Note=Other end matches $hn1"), "\n";
    print join("\t", $hn2, 'dundeex', 'BES', $hs2, $he2,      $sc2, '+', '.',
	       "ID=$gi2;Parent=$clone.TV;Note=TV"), "\n";
  }
  if($st2 eq 'C'){
    print join("\t", $hn2, 'dundeex', 'BAC', $hs2-1000, $he2,  '.', '+', '.',
	       "ID=$clone.TV;Name=$clone.TV;Note=Other end matches $hn1"), "\n";
    print join("\t", $hn2, 'dundeex', 'BES', $hs2, $he2,      $sc2, '-', '.',
	       "ID=$gi2;Parent=$clone.TV;Note=TV"), "\n";
  }
}

warn "OK\n";



sub hits_by_score {
  # score
  $b->[0] <=> $a->[0] ||
    # ident
    $b->[9] <=> $a->[9];
}

