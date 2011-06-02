#!/usr/bin/perl -w

use strict;

=head1 DESCRIPTION

This script produces GFF from specific 454 PE alignment summaries. The
script removes very short or very long pairs, so they will never be
seen!

=cut

use Getopt::Long;

## Set the source field, if desired
my $source_tag = 'dundee';

## Do we have access to trim file(s)?
my @trim_files;

## How about read file(s)?
my @read_files;

## Set to 1 to enable debugging output
my $verbose = 0;

GetOptions( "source=s"  => \$source_tag,
            "v"         => \$verbose,
            "trim=s{,}" => \@trim_files,
            "read=s{,}" => \@read_files,
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
    ## Actually all we care about is the length of the read
    $trim_data{$trim[0]} = $trim[2]
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
    $read_data{$read[0]} = [@read[1,2,3]];
  }
}

warn "found ", scalar keys %read_data, "\n"
  if $verbose > 0;





## Read in the alignment results

warn "\nparsing results\n";

while(<>){
  
  my ($clone, $status, $distance,
      $hn1, $hs1, $st1,
      $hn2, $hs2, $st2) = split/\s+/;
  
  # debugging
  die "failed to find data for $clone\n"
    unless
      exists $trim_data{"$clone\_left" } &&
      exists $trim_data{"$clone\_right"} &&
      exists $read_data{"$clone\_left" } &&
      exists $read_data{"$clone\_right"};
  
  
  
  ## Read end points are not recorded in the pair file
  my $he1 = $hs1 + ($trim_data{"$clone\_left" } || die ); # 400 );
  my $he2 = $hs2 + ($trim_data{"$clone\_right"} || die ); # 400 );
  
  ## Get additional read mapping data
  my ($ms1, $ma1, $pm1) = @{$read_data{"$clone\_left" }};
  my ($ms2, $ma2, $pm2) = @{$read_data{"$clone\_right"}};
  
  ## Only want 'Full' matches
  next if $ms1 ne "Full";
  next if $ms2 ne "Full";
  
  
  
  ## Process pairs
  
  if($status eq 'TruePair'){
    print
      join("\t",
           $clone, $status, $distance,
           $hn1, $hs1, $st1,
           $hn2, $hs2, $st2
          ), "\n"
            if $verbose > 1;
    
    ## Determine the direction of the clone on the sequence
    
    my $direction;
    
    if(   $st1 eq '+' && $st2 eq '-'){
      ## Sanity check
      die "what 1?\n" if $hs1 > $hs2;
      $direction = +1;
    }
    elsif($st1 eq '-' && $st2 eq '+'){
      ## Sanity check
      die "what 2?\n" if $hs1 < $hs2;
      $direction = -1;
    }
    else{
      ## Sanity check
      die "what 3? : $_\n";
    }
    
    ## One last sanity check
    if($hn1 ne $hn2){
      die;
    }
    
    
    
    ## Filter clones by distance:
    
    ## Check their distance
    print "I:\tdistance is $distance\n"
      if $verbose > 1;
    
    ## Give the clone the lowest 'mapping accuracy' of the mate-pair
    my $ma0 = (sort($ma1, $ma2))[0];
    
    if ($direction == +1){
      ## OK! Print three lines of GFF
      print join("\t", $hn1, $source_tag, 'clone',     $hs1, $he2, $ma0, '+', '.', "ID=$clone;Name=$clone;dist=$distance"), "\n";
      print join("\t", $hn1, $source_tag, 'clone end', $hs1, $he1, $ma1, '+', '.', "ID=$clone\_left;Parent=$clone"), "\n";
      print join("\t", $hn1, $source_tag, 'clone end', $hs2, $he2, $ma2, '-', '.', "ID=$clone\_right;Parent=$clone"), "\n";
    }
    if ($direction == -1){
      ## OK! Print three lines of GFF
      print join("\t", $hn1, $source_tag, 'clone',     $hs2, $he1, $ma0, '-', '.', "ID=$clone;Name=$clone;dist=$distance"), "\n";
      print join("\t", $hn1, $source_tag, 'clone end', $hs1, $he1, $ma1, '-', '.', "ID=$clone\_left;Parent=$clone"), "\n";
      print join("\t", $hn1, $source_tag, 'clone end', $hs2, $he2, $ma2, '+', '.', "ID=$clone\_right;Parent=$clone"), "\n";
    }
  }
  
  
  
  if($status eq 'FalsePair'){
    
    print
      join("\t",
           $clone, $status, $distance,
           $hn1, $hs1, $st1,
           $hn2, $hs2, $st2
          ), "\n"
            if $verbose > 1;
    
    ## Sanity checks
    if($hn1 eq $hn2){
      ## Orientation or distane error (not interested here)
      next;
    }
    
    
    
    ## Remeber 1 = left
    ## Remeber 2 = right

    ## ONE
    if(0){}
    elsif($st1 eq '+'){
      print join("\t", $hn1, "$source_tag link", 'clone',     $hs1, $he1+1000, $ma1, '+', '.', "ID=$clone\_left;Name=$clone\_left;Note=Other end matches $hn2"), "\n";
      print join("\t", $hn1, "$source_tag link", 'clone end', $hs1, $he1,      $ma1, '+', '.', "Parent=$clone\_left"), "\n";
    }
    elsif($st1 eq '-'){
      print join("\t", $hn1, "$source_tag link", 'clone',     $hs1-1000, $he1, $ma1, '-', '.', "ID=$clone\_left;Name=$clone\_left;Note=Other end matches $hn2"), "\n";
      print join("\t", $hn1, "$source_tag link", 'clone end', $hs1,      $he1, $ma1, '-', '.', "Parent=$clone\_left"), "\n";
    }
    else{die}
    
    ## TWO
    if(0){}
    elsif($st2 eq '+'){
      print join("\t", $hn2, "$source_tag link", 'clone',     $hs2, $he2+1000, $ma2, '-', '.', "ID=$clone\_right;Name=$clone\_right;Note=Other end matches $hn1"), "\n";
      print join("\t", $hn2, "$source_tag link", 'clone end', $hs2, $he2,      $ma2, '+', '.', "Parent=$clone\_right"), "\n";
    }
    elsif($st2 eq '-'){
      print join("\t", $hn2, "$source_tag link", 'clone',     $hs2-1000, $he2, $ma2, '+', '.', "ID=$clone\_right;Name=$clone\_right;Note=Other end matches $hn1"), "\n";
      print join("\t", $hn2, "$source_tag link", 'clone end', $hs2,      $he2, $ma2, '-', '.', "Parent=$clone\_right"), "\n";
    }
    else{die}
    
  }
  
}

warn "OK\n";
