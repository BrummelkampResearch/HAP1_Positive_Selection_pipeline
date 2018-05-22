#! /usr/local/bin/perl -w

# purpose:
# parsed files and create data for fisher test
# input files:
#     coupledtogenes_correct+symbol
#     coupledtogenes_incorrect+symbol
#          
# format: tab delimited 
# chr10   100195839       100195872       Choleratoxin    HPS1
#
# Output:
# ABT     AATF    1       492     2       31541
#
# author: Bingbing Yuan
#

use strict;

my ($HIGH_FILE, $LOW_FILE, $SCREENNAME) = @ARGV;

# only screan name and symbol columns are used, ignore the rest of column
my $SCREEN_COL = 3;
my $SYMBOL     = 4;

# 
# read correct_orentation file and save into hash of hash
# key_1=screenName
# for each screenName: key_2=symbol value=counts with this symbol

my ($HIGH_HREF, $HIGH_COUNT_HREF) =  fileToHashOfHash($HIGH_FILE);

# read the incorrect_orientation file (control), and save it into simple hash table 
# because screen name is not considered
# key=symbol; value=counts with this symbol

my ($LOW_HREF, $LOW_TOTAL) = fileToHash($LOW_FILE);


# 
# create the 4 values for fisher test
#
foreach my $screen (sort keys %$HIGH_HREF) {
  my $count_per_screen = $HIGH_COUNT_HREF->{$screen};

  foreach my $symbol (sort keys %{$HIGH_HREF->{$screen}} ) {
    my @data = ();
    my $other_symbol = $count_per_screen - $HIGH_HREF->{$screen}->{$symbol};

    my $control_symbol = 0;
    if ( $LOW_HREF->{$symbol} ) {
      $control_symbol = $LOW_HREF->{$symbol};
    }

    push (@data, $screen, $symbol, $HIGH_HREF->{$screen}->{$symbol}, $other_symbol, $control_symbol, $LOW_TOTAL-$control_symbol);
    print join("\t", @data), "\n";
  
#    print "$HIGH_HREF->{$screen}->{$symbol}, count_per_screen=$count_per_screen, $control_symbol, $LOW_TOTAL \n";
    
  }
}


#sub get_screen_count {
#  my $href = shift;
#  my $total = 0;

#  foreach my $key (keys %$href) {
#    print "symbol=$key\n";
#    $total = $total + $href->{$key};
#  }
#  return $total;
#}



sub fileToHash {
  my $file = shift;
  my %hash;
  my $total = 0;
  
  open(FIRSTFILE, $file) || die "can not read $file $!\n";
  while (<FIRSTFILE>) {
    chomp();
    next if ($_ =~ /^\#/);
    next if ($_ =~ /^\s*$/);
    my @arr = split(/\t/, $_);    
    
    my $_symbol      = $arr[$SYMBOL];
    
    if (! $hash{$_symbol}) {
      $hash{$_symbol} = 0;
    }
    $hash{$_symbol} ++;
    $total++;
  }
  close(FIRSTFILE);

  return (\%hash, $total); 
}



sub fileToHashOfHash {
  my $file = shift;
  my (%hash, %hash_count);
  
  open(FIRSTFILE, $file) || die "can not read $file $!\n";
  while (<FIRSTFILE>) {
    chomp();
    next if ($_ =~ /^\#/);
    next if ($_ =~ /^\s*$/);
    my @arr = split(/\t/, $_);    
    
    my $_screen_name = $arr[$SCREEN_COL];
    my $_symbol      = $arr[$SYMBOL];

    if ( ! $hash{$_screen_name}{$_symbol} ) {
      $hash{$_screen_name}{$_symbol} = 0;
    }
    $hash{$_screen_name}{$_symbol} ++;


    if (! $hash_count{$_screen_name} ) {
      $hash_count{$_screen_name} = 0;
    }
    $hash_count{$_screen_name} ++;

  }
  close(FIRSTFILE);

  return (\%hash, \%hash_count); 
}
