#!/usr/bin/perl -w

use strict;

=head1 DESCRIPTION

Thsi script produces GFF from BES alignments. It uses a specific
regexp to match one end to another, and doesn't check the BES
separation distance on the target sequence.

=cut


## We expect to be passed the SSAHA2 alignment file for the BES
die "pass ssaha file\n"
  unless @ARGV;



## Set to 1 to enable debugging output
my $verbose = 0;



## Read in the sequence data table so that we can match the sequence
## identifier in the alignment file to the more informative sequence
## header (parsed into the table elsewhere).

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



## Additionally we read in a list of BES that get filtered

warn "\nparsing filtered list\n";

my $filt_list_file = "filt_bes.list";

open F, '<', $filt_list_file
  or die "failed to open $filt_list_file : $! \n";

my %filt_list;

while(<F>){
  chomp;
  $filt_list{$_}++
}

warn "got ", scalar keys %filt_list, "\n";







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
  my $gi1 = $sequence_pairs{$clone}{'TP'} || 'seq missing';
  my $gi2 = $sequence_pairs{$clone}{'TV'} || 'seq missing';
  
  my @hits1 = exists($hits{$gi1}) ? @{$hits{$gi1}} : ();
  my @hits2 = exists($hits{$gi2}) ? @{$hits{$gi2}} : ();
  
  ## Debugging
  print "G:\t$gi1\t", scalar(@hits1), "\n"
    if $verbose > 0;
  print "G:\t$gi2\t", scalar(@hits2), "\n"
    if $verbose > 0;
  
  
  
  
  
  ## Create some failure flags?
  my ($missing_f, $rep_f, $ali_f) = (0, 0, 0);
  
  if ($gi1 eq 'seq missing' ||
      $gi2 eq 'seq missing'){
    $missing_f++;
  }
  
  if (@hits1 == 0){
    exists( $filt_list{ $sequence_pairs{$clone}{'TP'} || 'seq missing' } ) ?
      $rep_f++ : $ali_f++;
  }
  if (@hits2 == 0){
    exists( $filt_list{ $sequence_pairs{$clone}{'TV'} || 'seq missing' } ) ?
      $rep_f++ : $ali_f++;
  }
  
  ## OK
  
  
  
  
  
  ## Try to pair the hits (best hits first)
  
  ## NOTE: if the best hits don't form a 'good pair', we give up on
  ## the clone.

  ## NOTE: about 3% of the clones could be 'recovered' by hunting for
  ## a lower scoring 'good pair', however, their quality is dubious,
  ## and they are complex to handle! We skip them here.
  
  my ($scaff_f, $ori_f, $short_f, $long_f) = (0, 0, 0, 0);
  
 HIT:
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
	my $bes_distance = $he2 - $hs1 + 1;
	print "I:\tdistance is $bes_distance\n"
	  if $verbose > 0;
	
	if ($bes_distance < 5000){
	  ## Debugging
	  print "W:\treally incorrect orientation\n"
	    if $verbose > 0;
	  $short_f++;
	  last HIT;
	}
	
	if ($bes_distance > 400_000){
	  ## Debugging
	  print "W:\ttoo long!\n"
	    if $verbose > 0;
	  $long_f++;
	  last HIT;
	}
	
	## OK! Print three lines of GFF
	print join("\t", $hn1, 'dundee', 'BAC', $hs1, $he2,  '.', '+', '.', "ID=$clone;Name=$clone"), "\n";
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
	  $ori_f++;
	  last HIT;
	}
	
	## OK. Check their distance
	my $bes_distance = $he1 - $hs2 + 1;
	print "I:\tdistance is $bes_distance\n"
	  if $verbose > 0;
	
	if ($bes_distance < 1000){
	  ## Debugging
	  print "W:\treally incorrect orientation\n"
	    if $verbose > 0;
	  $short_f++;
	  last HIT;
	}
	
	if ($bes_distance > 400_000){
	  ## Debugging
	  print "W:\ttoo long!\n"
	    if $verbose > 0;
	  $long_f++;
	  last HIT;
	}
	
	## OK! Print three lines of GFF
	print join("\t", $hn1, 'dundee', 'BAC', $hs2, $he1,  '.', '-', '.', "ID=$clone;Name=$clone"), "\n";
	print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1, $sc1, '-', '.', "ID=$gi1;Parent=$clone;Note=TP"), "\n";
	print join("\t", $hn1, 'dundee', 'BES', $hs2, $he2, $sc2, '+', '.', "ID=$gi2;Parent=$clone;Note=TV"), "\n";
	
	## DONE!
	next CLONE;
      }
      
    }
  }
  
  
  
  
  
  
  
  ## Nothing paired, falling through to single end printer\n";
  
  ## NB! We simply don't care about a whole class of failed pairs!
  next if $missing_f || $rep_f || $ali_f;
  
  ## NB! We simply don't care about another class of failed pairs!
  next if $ori_f || $short_f || $long_f;
  
  
  
  ## The only remaining unpaired ends are those that span
  ## superscaffolds!
  
#   ## Remeber 1 = TP = forward
#   for my $hit1 (sort hits_by_score @hits1){
    
#     my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1,
# 	$al1, $nt1, $ql1) = @$hit1;
    
#     ## Debugging
#     print join("\t", '1:', @$hit1), "\n"
#       if $verbose > 0;
    
#     if($st1 eq 'F'){
#       print join("\t", $hn1, 'dundee', 'BAC', $hs1, $he1+1000,  '.', '+', '.', "ID=$clone.TP;Name=$clone.TP"), "\n";
#       print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1,      $sc1, '+', '.', "ID=$gi1;Parent=$clone.TP;Note=TP"), "\n";
#     }
#     if($st1 eq 'C'){
#       print join("\t", $hn1, 'dundee', 'BAC', $hs1-1000, $he1,  '.', '-', '.', "ID=$clone.TP;Name=$clone.TP"), "\n";
#       print join("\t", $hn1, 'dundee', 'BES', $hs1, $he1,      $sc1, '-', '.', "ID=$gi1;Parent=$clone.TP;Note=TP"), "\n";
#     }
#     last;
#   }
  
#   ## Remeber 2 = TV = reverse
#   for my $hit2 (sort hits_by_score @hits2){
    
#     my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2,
# 	$al2, $nt2, $ql2) = @$hit2;
    
#     ## Debugging
#     print join("\t", '2:', @$hit2), "\n"
#       if $verbose > 0;
    
#     if($st2 eq 'F'){
#       print join("\t", $hn2, 'dundee', 'BAC', $hs2, $he2+1000,  '.', '-', '.', "ID=$clone.TV;Name=$clone.TV"), "\n";
#       print join("\t", $hn2, 'dundee', 'BES', $hs2, $he2,      $sc2, '+', '.', "ID=$gi2;Parent=$clone.TV;Note=TV"), "\n";
#     }
#     if($st2 eq 'C'){
#       print join("\t", $hn2, 'dundee', 'BAC', $hs2-1000, $he2,  '.', '+', '.', "ID=$clone.TV;Name=$clone.TV"), "\n";
#       print join("\t", $hn2, 'dundee', 'BES', $hs2, $he2,      $sc2, '-', '.', "ID=$gi2;Parent=$clone.TV;Note=TV"), "\n";
#     }
#     last;
#   }
  
  
  
  ## Remeber 1 = TP = forward
  ## Remeber 2 = TV = reverse
  my $hit1 = (sort hits_by_score @hits1)[0];
  my $hit2 = (sort hits_by_score @hits2)[0];
  
  my ($sc1, $qn1, $hn1,  $qs1, $qe1, $hs1, $he1, $st1, $al1, $nt1, $ql1) = @$hit1;
  my ($sc2, $qn2, $hn2,  $qs2, $qe2, $hs2, $he2, $st2, $al2, $nt2, $ql2) = @$hit2;
  
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

