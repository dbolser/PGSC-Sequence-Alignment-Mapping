#!/usr/bin/perl -w

use strict;

die "pass ssaha file\n"
  unless @ARGV;



## Set to 1 to enable debugging output
my $verbose = 0;



## Parse sequence data table

warn "\nparsing sequence data table\n";

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
  
  ## No direction shoudl be duplicated!
  die if
    exists($sequence_pairs{$clone}{$direction});
  
  ## Store the 'putative pair'
  $sequence_pairs{$clone}{$direction} = $gi;
}

warn "got $seq_count BES (and ",
  scalar keys %sequence_pairs, " BACs)\n";

warn "\nNote, not BACs have BES for both directions!\n";



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





## Produce the GFF, paired and all

warn "\ndumping GFF\n";

my %paired;

CLONE:foreach my $clone (keys %sequence_pairs){
  
  ## Debugging
  print "C:\t$clone\n"
    if $verbose > 0;
  
  ## Get the hits for both sequences
  my $gi1 = $sequence_pairs{$clone}{'TP'} || 'missing';
  my $gi2 = $sequence_pairs{$clone}{'TV'} || 'missing';
  
  my @hits1 = exists($hits{$gi1}) ? @{$hits{$gi1}} : ();
  my @hits2 = exists($hits{$gi2}) ? @{$hits{$gi2}} : ();
  
  ## Debugging
  print "G:\t$gi1\t", scalar(@hits1), "\n"
    if $verbose > 0;
  print "G:\t$gi2\t", scalar(@hits2), "\n"
    if $verbose > 0;
  
  
  
  ## Try to pair the hits (best hits first)
  
  for my $hit1 (sort hits_by_score @hits1){
    for my $hit2 (sort hits_by_score @hits2){
      
      my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1,
	  $al1, $nt1, $ql1) = @$hit1;
      my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2,
	  $al2, $nt2, $ql2) = @$hit2;
      
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
	next;
      }
      
      ## Check orientation 1/2
      if ($st1 eq 'F'){
	if ($st2 ne 'C'){
	  ## Debugging
	  print "W:\tincorrect orientation\n"
	    if $verbose > 0;
	  next;
	}
	
	## OK, print two lines of GFF
	my $bes_distance = $he2 - $hs1 + 1;
	
	if ($bes_distance < 1000){
	  ## Debugging
	  print "W:\treally incorrect orientation\n"
	    if $verbose > 0;
	  next;
	}
	
	print join("\t", $hn1, 'dundee', 'BAC', $hs1, $he2,  '.', '+', '.', "ID=$clone"), "\n";
	print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1, $sc1, '+', '.', "ID=$gi1;Parent=$clone;Note=TP"), "\n";
	print join("\t", $hn1, 'dundee', 'BES', $hs2, $he2, $sc2, '-', '.', "ID=$gi2;Parent=$clone;Note=TV"), "\n";
	
	## DONE!
	next CLONE;
      }
      
      ## Check orientation 2/2
      if ($st1 eq 'C'){
	if ($st2 ne 'F'){
	  ## Debugging
	  print "W:\tincorrect orientation\n"
	    if $verbose > 0;
	  next;
	}
	
	## OK, print two lines of GFF
	my $bes_distance = $he1 - $hs2 + 1;
	
	if ($bes_distance < 1000){
	  ## Debugging
	  print "W:\treally incorrect orientation\n"
	    if $verbose > 0;
	  next;
	}
	
	print join("\t", $hn1, 'dundee', 'BAC', $hs2, $he1,  '.', '-', '.', "ID=$clone"), "\n";
	print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1, $sc1, '-', '.', "ID=$gi1;Parent=$clone;Note=TP"), "\n";
	print join("\t", $hn1, 'dundee', 'BES', $hs2, $he2, $sc2, '+', '.', "ID=$gi2;Parent=$clone;Note=TV"), "\n";
	
	## DONE!
	next CLONE;
      }
    }
  }
  
  
  
  ## Nothing paired, falling through to single end printer\n";
  
  ## Remeber 1 = TP = forward
  for my $hit1 (sort hits_by_score @hits1){
    
    my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1,
	$al1, $nt1, $ql1) = @$hit1;
    
    ## Debugging
    print join("\t", '1:', @$hit1), "\n"
      if $verbose > 0;
    
    if($st1 eq 'F'){
      print join("\t", $hn1, 'dundee', 'BAC', $hs1, $he1+1000,  '.', '+', '.', "ID=$clone.TP"), "\n";
      print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1,      $sc1, '+', '.', "ID=$gi1;Parent=$clone.TP;Note=unpaired TP"), "\n";
    }
    if($st1 eq 'C'){
      print join("\t", $hn1, 'dundee', 'BAC', $hs1-1000, $he1,  '.', '-', '.', "ID=$clone.TP"), "\n";
      print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1,      $sc1, '-', '.', "ID=$gi1;Parent=$clone.TP;Note=unpaired TP"), "\n";
    }
    last;
  }
  
  ## Remeber 2 = TV = reverse
  for my $hit2 (sort hits_by_score @hits2){
    
    my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2,
	$al2, $nt2, $ql2) = @$hit2;
    
    ## Debugging
    print join("\t", '2:', @$hit2), "\n"
      if $verbose > 0;
    
    if($st2 eq 'F'){
      print join("\t", $hn2, 'dundee', 'BAC', $hs2, $he2+1000,  '.', '-', '.', "ID=$clone.TV"), "\n";
      print join("\t", $hn2, 'dundee', 'BES', $hs2, $he2,      $sc2, '+', '.', "ID=$gi2;Parent=$clone.TV;Note=unpaired TV"), "\n";
    }
    if($st2 eq 'C'){
      print join("\t", $hn2, 'dundee', 'BAC', $hs2-1000, $he2,  '.', '+', '.', "ID=$clone.TV"), "\n";
      print join("\t", $hn2, 'dundee', 'BES', $hs2, $he2,      $sc2, '-', '.', "ID=$gi2;Parent=$clone.TV;Note=unpaired TV"), "\n";
    }
    last;
  }
}

warn "OK\n";



sub hits_by_score {
  # score
  $b->[0] <=> $a->[0] ||
    # ident
    $b->[9] <=> $a->[9];
}
