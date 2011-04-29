#!/usr/bin/perl -w

use strict;

use Data::Dumper;

use Getopt::Long;


## Mates file for PEMP data
my $mates_file;

## Set to 1 to enable debugging output
my $verbose = 0;



GetOptions(
           'mates_file=s'   => \$mates_file,
           'verbose+'       => \$verbose,
          )
  or die "failed to communicate!\n";

## SANITY CHECK COMMAND LINE OPTIONS

die "please pass a mates file\n"
  unless defined($mates_file) && -s $mates_file;

## We expect to be passed some GFF

die "pass GFF\n"
  unless @ARGV;



## Process the mates file...
my (%mates, $curr_lib);

open M, '<', $mates_file
  or die "failed to open file '$mates_file' : $!\n";

while(<M>){
  chomp;
  if(/^library/){
    my (undef, $lib, $min, $max, $lib_regexp) = split/\t/;
    $curr_lib = $lib;
    
    $mates{$curr_lib}{min} = $min;
    $mates{$curr_lib}{max} = $max;
    $mates{$curr_lib}{lib_regexp} = qr/$lib_regexp/;
  }
  if(/^pair/){
    my (undef, $forward_regexp, $reverse_regexp) = split/\t/;
    
    $mates{$curr_lib}{forward_regexp} = qr/$forward_regexp/;
    $mates{$curr_lib}{reverse_regexp} = qr/$reverse_regexp/;
  }
}

warn "Read ", scalar keys %mates, " mates libs\n";
warn "OK\n\n";





## Read in GFF, pair ends and match library

warn "processing GFF\n";

my (%hits, @gff_directives);

while(<>){
  chomp;
  
  ## Process GFF directives
  if (/^#/){
    if(/^##/){
      push @gff_directives, $_;
    }
    next;
  }
  
  
  ## Process GFF
  my ($seq_id, $source, $type, $start, $end, $score,
      $strand, $phase, $args
     ) = split/\t/;
  
  ## Process 'args'
  my %args;
  for (split(/;/, $args)){
    my ($key, $value) = split/=/;
    die if exists $args{$key};
    $args{$key} = $value;
  }
  
  my $id = $args{'ID'};
  
  
  
  ## Match this 'clone end' id to a library
  my @libs =
    grep {$id =~ $mates{$_}{lib_regexp}} keys %mates;
  
  ## One library only please!
  die "failed to match a unique library for '$id'\nmatched: '",
    join("'\t'", @libs), "'\n" if @libs != 1;
  
  
  
  ## Prolly safer to do this with the library regexp!
  my $base_id = substr($id, 0, -2);
  my $direction = substr($id, -1);
  
  push @{$hits{$libs[0]}{$base_id}{$direction}},
    [$seq_id, $source, $type, $start, $end, $score,
     $strand, $phase, \%args];
}

## Make more informative?
warn "got ", scalar keys %hits, "\n";





## Produce the GFF

warn "\ndumping GFF\n";

print "$_\n" for @gff_directives;



my %failure_code;

foreach my $lib (keys %hits){
  print "L:\t$lib\n"
    if $verbose > 0;
  
  foreach my $clone (keys %{$hits{$lib}}){
    ## Debugging
    print "C:\t\t$clone\n"
      if $verbose > 1;
    
    ## Get the hits for both sequences
    my $hits1 = $hits{$lib}{$clone}{'F'} || []; # FORWARD
    my $hits2 = $hits{$lib}{$clone}{'R'} || []; # REVERSE
    
    ## Debugging
    print "G:\t\t$clone\_F\t", scalar(@$hits1), "\n"
      if $verbose > 1;
    print "G:\t\t$clone\_R\t", scalar(@$hits2), "\n"
      if $verbose > 1;
    
    
    
    ## Failure accounting
    
    ## Hits for both ends? If not, we give up now (we don't want to
    ## 'pollute' the GFF with these kinds of problem, as they can be
    ## misleading... I guess).
    
    if   ( @$hits1 == 0 && @$hits2 == 0 ){
      $failure_code{'no align'}++;
      next;
    }
    elsif( @$hits1 == 0 || @$hits2 == 0 ){
      $failure_code{'end no align'}++;
      next;
    }
    
    
    
    ## Try to pair the hits (best hits first)
    
    ## NOTE: if the best hits for each end don't form a 'good pair',
    ## we give up on the clone. About 3% of the clones could be
    ## 'recovered' by hunting for a lower scoring 'good pair',
    ## however, their quality is dubious, and they are complex to
    ## handle! We skip them here.
    
    my $hit1 = (sort hits_by_score @$hits1)[0]; # FORWARD
    my $hit2 = (sort hits_by_score @$hits2)[0]; # REVERSE
    
    my ($seq_id1, $source1, $type1, $start1, $end1, $score1,
        $strand1, $phase1, $args_ref1) = @$hit1;
    my ($seq_id2, $source2, $type2, $start2, $end2, $score2,
        $strand2, $phase2, $args_ref2) = @$hit2;
    
    ## Debugging
    print join("\t", '1:', @$hit1), "\n"
      if $verbose > 2;
    print join("\t", '2:', @$hit2), "\n"
      if $verbose > 2;
    
    
    
    ## Check they hit the same scaffold. If not they must span a pair
    ## of scaffolds, and are handled differently (below).
    
    ##
    ## PAIRED!
    ##
    
    if ($seq_id1 eq $seq_id2){
      
      ## Calculate clone size, sidestepping issues of incorrect or
      ## incorrect orientation
      my ($st, $en) =
        range($start1, $end1,
              $start2, $end2);
      my $insert_size = $en - $st;
      
      ## Check orientation and establish clone orientation
      my $good_pair = 'b';  ## b = bad
      my $clone_strand = 0; ## 0 = unknown
      
      if($start1 < $start2){
        if($strand1 eq '+' && $strand2 eq '-'){
          $good_pair = 'g'; ## g = good
          $clone_strand = '+'; ## a +ve clone
        }
      }
      else{
        if($strand2 eq '+' && $strand1 eq '-'){
          $good_pair = 'g'; ## g = good
          $clone_strand = '-'; ## a -ve clone
        }
      }
      ## a 'bad' pair, clone strand is undefined
      
      ## Check length
      my $happy;
      $happy = 'h';
      $happy = 's' if $insert_size < $mates{$lib}{min};     ## short
      $happy = 'l' if $insert_size > $mates{$lib}{max};     ## long
      $happy = 'f' if $insert_size > $mates{$lib}{max} * 4; ## fail
      
      print
        join("\t",
             $seq_id1, "$source1 $happy$good_pair", 'clone', $st, $en, '.',
             $clone_strand, '.',
             join(';',
                  "ID=$clone", "Name=$clone",
                 ),
            ), "\n";
      
      print
        join("\t",
             $seq_id1, $source1, $type1, $start1, $end1, $score1,
             $strand1, $phase1,
             join(';',
                  "Parent=$clone", map("$_=". $args_ref1->{$_}, keys %$args_ref1)
                 ),
            ), "\n";
      
      print
        join("\t",
             $seq_id2, $source2, $type2, $start2, $end2, $score2,
             $strand2, $phase2,
             join(';',
                  "Parent=$clone", map("$_=". $args_ref2->{$_}, keys %$args_ref2)
                 ),
            ), "\n";
      
      ## Failure accounting
      $failure_code{"paired: $source1 $happy$good_pair"}++;
    }
    
    
    
    ##
    ## DID NOT PAIR: PRINT EACH CLONE END SEPARATELY
    ##
    
    else{
      print
        join("\t",
             $seq_id1, "$source1 l", 'clone',
             $start1 + ($strand1 eq '+' ? 0 : -10000),
             $end1   + ($strand1 eq '+' ? +10000 : 0),
             $score1, $strand1, '.',
             join(';',
                  "ID=$clone\_F", "Name=$clone\_F",
                  "Note=Other end matches ". $seq_id2,
                 ),
            ), "\n";
      
      ## ID tweek to avoid ID overlap
      $args_ref1->{'ID'} = '_'. $args_ref1->{'ID'};
      
      print
        join("\t",
             $seq_id1, "$source1", $type1, $start1, $end1, $score1,
             $strand1, '.',
             join(';',
                  "Parent=$clone\_F", map("$_=". $args_ref1->{$_}, keys %$args_ref1)
                 ),
            ), "\n";
      
      print
        join("\t",
             $seq_id2, "$source2 l", 'clone',
             $start2 + ($strand2 eq '+' ? 0 : -10000),
             $end2   + ($strand2 eq '+' ? +10000 : 0),
             $score2, $strand2, '.',
             join(';',
                  "ID=$clone\_R", "Name=$clone\_R",
                  "Note=Other end matches ". $seq_id1,
                 ),
            ), "\n";
      
      ## ID tweek to avoid ID overlap
      $args_ref2->{'ID'} = '_'. $args_ref2->{'ID'};
      
      print
        join("\t",
             $seq_id2, "$source2", $type2, $start2, $end2, $score2,
             $strand2, '.',
             join(';',
                  "Parent=$clone\_R", map("$_=". $args_ref2->{$_}, keys %$args_ref2)
                 ),
            ), "\n";
      
    }
  }
}

## output failures
warn "$_\t$failure_code{$_}\n" for sort keys %failure_code;

warn "OK\n";





sub hits_by_score {
  # score
  $b->[5] <=> $a->[5] ||
    # target
    $b->[0] cmp $a->[0];
}



sub range {
  my $min = $_[+0];
  my $max = $_[-1];
  
  for(@_[1..($#_-1)]){
    $min = $_ if $_ < $min;
    $max = $_ if $_ > $max;
  }
  
  return ($min, $max);
}
