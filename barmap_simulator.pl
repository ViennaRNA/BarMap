#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2022-09-21 21:17:09 ivo>

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use strict;
use warnings;
use Carp;
use List::Util qw(all);
use Scalar::Util qw(reftype looks_like_number);
use YAML;
use vars qw/
  $T0 $T8 $TX $P0 $TINC $SIM_TIMES_FILE $TEMP $EQ @FILES $TREEKIN %ENV
  $RATESUFFIX $RECOVER
/;

# defaults for global(s)
$TREEKIN = 'treekin';
$T0 = 0.1;
$T8 = 10;
$TX = -1;
$TINC = 1.02;
$P0 = "";
$EQ = 0;
$TEMP = 37.0;
$RATESUFFIX="rates.bin";
my $outfile=undef;
my $OUT;

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1)
    unless GetOptions(
          "t0=f"  => \$T0,
		      "t8=f"  => \$T8,
		      "tX=f"  => \$TX,
		      "inc=f" => \$TINC,
		      "times-file=s" => \$SIM_TIMES_FILE,
		      "eq=i"  => \$EQ,
		      "T=f"   => \$TEMP,
		      "p0=s"  => \&process_opt_p0,
		      "rs=s"  => \$RATESUFFIX,
		      "o=s"   => \$outfile,
		      "man"   => sub{pod2usage(-verbose => 2)},
		      "help"  => sub{pod2usage(-verbose => 1)},
          "recover" => \$RECOVER,
       );

defined $outfile ? open($OUT, ">", $outfile) : open($OUT, ">&STDOUT");
die $! unless (defined $OUT);

my ($rs, $cs)  = parse_bar_map();
my $p_zeros    = $P0;
my $start_time = $T0;
my $stop_time  = $T8; 
my $sim_times_ref
  = defined $SIM_TIMES_FILE ? read_sim_times_file($SIM_TIMES_FILE) : {};
my $global_time = 0.0;
#print Dumper(\@FILES);

# main loop
my $flag = 0;
for (my $file = 0; $file<=$#FILES; $file++) {
  $stop_time = $TX, $flag = 1 if ($#FILES == $file) && ($TX > 0);
  print $OUT "# Rates file: ".$FILES[$file]." (#$file)\n";
  if ($EQ > 0){
    if ($file == $#FILES){
      $stop_time = $EQ;
      print $OUT "# Info: last landscape, setting t8 to $EQ\n";
    }
  }
#  print $OUT "p_zeros $p_zeros\n";
  my $command = build_command($p_zeros,
			      $stop_time,
			      $FILES[$file],
			      $sim_times_ref,
            $RECOVER,
            $RATESUFFIX,);
  print $OUT "# Cmd: $command\n";
  delete $sim_times_ref->{$FILES[$file]};       # delete used time params
  my ($stop, $densities, $tc) = do_simulation($command, $flag);
  if ($stop == -1) {
    print STDERR "Treekin ERROR at inputfile $FILES[$file]\n";
    last;
  }
  # dump time course to STDOUT
  dump_time_course($stop, $global_time, $tc, $flag);
  $global_time += $stop;

  last if $file == $#FILES;
  #print ">> file is $file<<\n";
  $p_zeros = remap_densities($densities, $cs->[$file], $cs->[$file+1]);
}
unlink qw(evecs.bin evals.bin) if $RECOVER;         # clean up

# Check if all simulation time parameters from the times file have been used.
my @unused_times_files = sort keys %$sim_times_ref;
warn 'WARNING: Encountered simulation time parameters for the following ',
     'unused rates files:', join ', ', @unused_times_files
  if @unused_times_files;

#---
sub dump_time_course {
  my ($stop, $global_time, $tc, $flag) = @_;

  # print time course
  foreach my $line (@$tc) {
    if ($line =~ m/^\#/) {
      print $OUT $line;
      next;
    }
    $line =~ s/^\s+//;
    my @values = split /\s+/, $line;
    my $time = shift @values;
    $time += $global_time;
    printf ($OUT "%22.20e ",$time);
    print $OUT join(" ", @values), "\n";
  }
  print $OUT "&\n";

  return if $flag;
  $stop += $global_time;
}

#---
sub calculate_stop_time {
  my $start = shift;
  my $time;
  for ($time=$start; $time <= $start+$T8; $time *= $TINC) {;}
  return $time;
}

#---
sub remap_densities {
  my ($densities, $current, $next) = @_;
  my $tuples = [ grep {!     (($_->[0] == -1)
			      || ($_->[1] == -1))} @{zip($current, $next)}];
  # there can be duplicate entries in the barmap file
  $tuples = uniq_tuples($tuples);
  #print STDERR ">>tuples\n";
  #print STDERR Dumper(\$tuples);
  #print STDERR ">>tuples DONE\n";
  my $dens_sum=0.;
  my @d = @$densities;
  shift @d;
  $dens_sum += $_ for @d;
  print STDERR "# we have ".eval($#$densities)." densities with a total density of $dens_sum from the previous landscape\n";

  # now get # of states in new landscape
  my $maxnew=0;
  for(my $i=0; $i<= $#$tuples; $i++) {
    my $new = $tuples->[$i]->[1];
    if($new > $maxnew){ $maxnew = $new }
  }

  my @p = (0) x ($maxnew);
  $maxnew=0;
  for(my $i=0; $i<= $#$tuples; $i++) {
    my $new = $tuples->[$i]->[1];
    my $old = $tuples->[$i]->[0];
    if ($old <= 0){
      print STDERR "skipping lmin $old\n";
      next;  # do we need a next at the end here?
    }
    $p[$new] += $densities->[$old];
    if($new > $maxnew){ $maxnew = $new }
  }
  my $total = 0.;
  $total += $_ for @p;

  print STDERR "# remapped! New total density is $total\n";
  print STDERR "# max lmin in new landscape is $maxnew\n\n";

  if (abs($total-1) > 1e-6){ # delta old/new population
    local $, = " ";
    print "\n";
    printf "%3d %3d  %7.5f %7.5f\n", @$_, $densities->[$_->[0]], $p[$_->[1]] for @$tuples;
    print "\n";
    croak "ERROR: Total remapped density is too large";
  }

  #print STDERR "\@p:\n";
  #print STDERR Dumper(\@p);

  my $yy = zip(number_sequence(), \@p);
  my $p0 = "";
  my $p0_sum = 0.0;
  for(my $i=1; $i<= $#p; $i++) {
    my $pop = $p[$i];
    next if ($pop == 0.);
    $p0 .= " --p0 $i=$pop";
      $p0_sum += $pop;
    }

  if ($p0_sum > 1.01) {
    print STDERR "$p0\np0_sum = $p0_sum\n";
    print STDERR Dumper([grep {!  (($_->[0] == -1)
			        || ($_->[1] == -1))} @$tuples]);
    print STDERR Dumper(\$p0);
    die $!
  }

  if ($p0_sum < 0.99) {
    print STDERR "$p0\np0_sum = $p0_sum\n";
    print STDERR Dumper([grep {!  (($_->[0] == -1)
			        || ($_->[1] == -1))} @$tuples]);
    print STDERR Dumper(\$p0);
    die $!
  }

  return $p0;
}

#---
sub uniq_tuples {
  my $tuples = shift;
  my @uniq = ();
  my $h = {map { $_->[0] => $_->[1] } @$tuples};
  foreach my $key (sort {$a <=> $b} keys %$h) {
    push @uniq, [$key, $h->{$key}];
  }
  return \@uniq;
}

#---
sub make_eigenval_links {
  # Recovering from pre-computed eigenvalues. Since Treekin reads these from
  # hardcoded file names 'evecs.bin' and 'evals.bin', links need to be created
  # from the respective eigenvalue file of the current rate matrix to match
  # the hardcoded file names.
  my ($rates_file, $rate_suffix) = @_;
  my $match_suffix = quotemeta $rate_suffix;    # escape regex meta chars
  my $base_name    = $rates_file =~ s/[.]$match_suffix//r;  # remove suffix
  my $evecs_file   = "$base_name.evecs.bin";    # eigenvecs / values are
  my $evals_file   = "$base_name.evals.bin";    # stored along w/ rate mat
  die "Eigenvector file '$evecs_file' and eigenvalue file '$evals_file' ",
      'required when using --recover, but not present or empty'
    unless -s $evecs_file and -s $evals_file;

  # Create symlinks.
  unlink qw(evecs.bin evals.bin);               # delete if existing
  symlink $evecs_file, 'evecs.bin';
  symlink $evals_file, 'evals.bin';
}

#---
sub build_command {
  my ($p_zeros, $t8_default, $rates_file, $times_ref, $recover, $ratesuffix) = @_;
  my $command = "$TREEKIN -m I";
  $command .= " --tinc=" . ($times_ref->{$rates_file}{inc  } // $TINC      );
  $command .= " --t0="   . ($times_ref->{$rates_file}{start} // $T0        );
  $command .= " --t8="   . ($times_ref->{$rates_file}{end  } // $t8_default);
  $command .= " --Temp=$TEMP";
  $command .= " --bin";
  if ($recover) {
    $command .= " --recoverE";
    make_eigenval_links $rates_file, $ratesuffix;
  }
  $command .= $p_zeros;
  $command .= " < $rates_file";
 # print Dumper($command);
  return $command;
}

#---
sub do_simulation {
  my ($command, $flag) = @_;
  print STDERR "$command\n\n";
  my @output = `$command`;

  # capture treekin error
  if ($#output == $[-1) {
    return (-1, undef, undef);
  }

  # get last line of time course not beginning with a '#' (i.e. no comment)
  my $lastline_index = (grep {$output[$_] !~ /^#/} 0..$#output)[-1];
  die 'No data found in Treekin output' unless defined $lastline_index;
  $output[$lastline_index] = '#' . $output[$lastline_index] if $flag == 1;
  my $lastline = $output[$lastline_index];

  $lastline =~ s/^\s+//;
  # extract values from last line of time course
  my @values = split /\s+/, $lastline;

  # check that we found at least two and nothing but numbers.
  die 'Received invalid last output line from Treekin: ', $lastline
    unless @values > 1 and all {looks_like_number $_} @values;

  my @valuescp = @values;
  shift @valuescp;
  my $vsum=0.;
  $vsum += $_ for @valuescp;
  print STDERR "total population at the end of the simulation is $vsum\n";
  # get stop time
  my $stop = $values[0];
  return ($stop, \@values, \@output);
}

#---
sub filter { [map { $_[0]->($_) ? $_ : ()} @{$_[1]}] }

#---
#sub is_defined { return 1 if defined $_->[1]; return 0 }

#---
sub process_opt_p0 {
  my $string = $_[1];
  if ($string =~ m/:/) {
    my @p0s = split /:/, $string;
    for my $p0 (@p0s) {
      my ($lmin, $value) = split(/=/, $p0);
      $P0 .= " --p0 $lmin=$value";
    }
  }
  else {
    my ($lmin, $value) = split(/=/, $string);
    $P0 .= " --p0 $lmin=$value";
  }
}

#---
sub parse_bar_map {
  my (@rows, $columns,$fn);
  while (<>) {
    s/^\#//, @FILES = split, next if $. == 1;
    chomp;
    # smash line into substrings of length 3 and capture only the odd ones
    # trim whitespace and convert empty strings to undef values
    # mtw 20171031: read up to 4 digits in group 1, ie allow max 9999 lmins
    # this breaks backward compatibility !!!
    my $row = [map {s/ //g; $_ eq "" ? -1 : $_ } (m/(....)(?:...)?/g)];
    # memorize row in a lol
    push @rows , $row;
    $columns = to_columns($#$row) if $. == 2;
    # memorize columns in a lol
    $columns->(zip(number_sequence(), $row));
  }
  foreach (@FILES){
    if ($RATESUFFIX eq "rates.r.bin"){
      $fn = fileparse($_, qr/\.r\.bar/);
    }
    else{
      $fn = fileparse($_, qr/\.bar/);
    }
    $_ = "$fn.".$RATESUFFIX;
    }
  # return lol's of rows and columns
  # print Dumper(@FILES);die;
  # print Dumper(\@rows);
  #print Dumper($columns->());
  return (\@rows, $columns->());
}

#---
# zip two lists togather into a list of tuples
# [a1, a2, ...], [b1, b2, ...] -> [[a1, b1], [a2, b2], ...]
sub zip {
  my ($A, $B) = @_;

  $A = array_to_iterator($A) if ref($A) eq 'ARRAY';
  $B = array_to_iterator($B) if ref($B) eq 'ARRAY';
 # print STDERR ">>zip() A:\n";
 # print STDERR Dumper($A->());
 # print STDERR ">>zip(): A DONE\n";
  my @result;
  while (   defined(my $a = $A->())
	 && defined(my $b = $B->())) {
    push @result, [$a, $b];
  }
 # print STDERR ">>zip() result:\n";
 # print STDERR Dumper(\@result);
 # print STDERR ">>zip(): result DONE\n";
  return \@result;
}

#---
# unzip a list of tupels into a tuple of lists
# [[a1, b1], [a2, b2], ...] -> [[a1, a2, ..], [b1, b2, ...]]
sub unzip {
  my (@A, @B);
  push(@A, $_->[0]), push(@B, $_->[1]) for @{shift()};

  return [\@A, \@B];
}

#---
# (val1, val2, ...) | [val1, val2, ...] -> {}
# converts an array into an iterator
sub array_to_iterator {
  my @array = (ref $_[0] eq 'ARRAY') ? @{shift()} : @_;

  return sub { shift @array }
}

#---
# n -> { }
# lazy version of infinit list [n..] of integers starting at value n
# implemented as closure
sub number_sequence {
  my $value = defined $_[0] ? $_[0] : 0;

  return sub { $value++ }
}

#---
# incrementally build lol from list of duples, implemented as closure
# [[idx1,val1], [idx2,val2], ...] -> [[values, val1],[values, val2], ...]
sub to_columns {
  my $lol = [ map{[]} $[..$_[0] ]; # empty lol with right dimension
  return sub {
    return $lol unless defined $_[0]; # return result
    map {push @{$lol->[$_->[0]]}, $_->[1]} @{$_[0]};
  }
}

#---
# Parses a YAML file containing individual simulation time parameters (i.e.
# start time, end time, and time increment) for one or more rates files.
# The file must contain a single dictionary where the keys are rates file
# names. The values are, again, single dictionaries which may contain one or
# more of the keys 'start', 'end', and 'inc', setting the respective parameter
# for that rates file.
#
# Returns a hash ref containing the described data structure.
sub read_sim_times_file {
  my ($times_yaml_file) = @_;

  # Parse simulation times from YAML file, skip undef entries (separators etc)
  my ($times_ref, @garbage)
    = grep {defined} YAML::LoadFile($times_yaml_file);

  # Check overall consistency.
  die "Simulation times file '$times_yaml_file' does not contain a ",
      "single hash ref"
    unless @garbage == 0
           and reftype $times_ref
           and reftype $times_ref eq reftype {}
           ;

  # Check consistency for each bar file.
  for my $rates_file (keys %$times_ref) {
    my $rates_times_ref = $times_ref->{$rates_file};
    die "Simulation times file '$times_yaml_file' does not contain ",
        "a hash ref for rates file '$rates_file'"
      unless reftype $rates_times_ref
             and reftype $rates_times_ref eq reftype {};

    # Check general validity of start, end, inc params and its values.
    my %is_valid_param = map {$_ => 1} qw(start end inc);
    while (my ($param, $value) = each %$rates_times_ref) {
      die "Simulation times file '$times_yaml_file' contains ",
          "invalid parameter '$param' for rates file '$rates_file'"
        unless $is_valid_param{$param};

      die "Simulation times file '$times_yaml_file' contains non-positive",
          " value '$value' for parameter '$param' of rates file '$rates_file'"
        unless defined $value and looks_like_number $value and $value > 0.;
    }
  }

  return $times_ref;
}


=head1 NAME

barmap_simulator.pl - perform a dynamic landscape simulation, e.g. for
co-transcriptional RNA folding

=head1 SYNOPSIS

barmap_simulator.pl [[-t0 I<FLOAT>] [-t8 I<FLOAT>] [-inc I<FLOAT>]
[-p0 I<STRING>] [-eq I<INT>] barmap-file

=head1 DESCRIPTION

This tool performs a BarMap simulation on a set of consecutive
landscapes. It expects a BarMap file that contains the exact landscape
mapping (pre-computed with I<bar_map.pl>) as last argument. Starting
from an initial population distribution on the first landscape (option
B<-p0>), this tool calls I<treekin> to simulate the folding dynamics
on each landscape for a certain time (measured in internal time,
option B<-t8>). At the end of such an increment, population of each
macrostate is mapped to the corrsponding macrostate in the next
landscape, and I<treekin> is called with this new start population.

A consistency check for overall popultion density is performed at
every elongation step. This ensures that no population density is lost
in the course of the simulation.

=head1 NOTES

I<treekin> simulations are performed on binary rates files per
default. This implies that for each C<.bar> file referenced in
C<barmap-file>, the corresponding binary rates file is present in the
current working directory.

=head1 OPTIONS

=over 4

=item B<-t0> I<FLOAT>

Set start time of kinetic simulation to I<FLOAT> (default: 0.1).

=item B<-t8> I<FLOAT>

Set stop time of kinetic simulation to I<FLOAT> (default: 10.0).

=item B<-inc> I<FLOAT>

Set time increment of kinetic simulation to I<FLOAT> (default: 1.02).

=item B<-times-file> I<STRING>

Path to a YAML file containing individual simulation time parameters (i.e.,
start time, end time, and time increment) for one or more rates files.
This is useful, e.g., to model a variable transcription rate along the RNA.
The values from the file overwrite the global values set via other arguments
(C<-t0>, C<-t8>, and C<-inc>) for the specified rates files.

The file is structured like this:

    "10.rates.bin":
        start: 1
        end: 100
        inc: 1.9
    "12.rates.bin":
        inc: 1.8

The order of the rates files is irrelevant.

=item B<-T> I<FLOAT>

Set simulation temperature to I<FLOAT> (default: 37.0).

=item B<-p0> I<STRING>

Specify initial population for minima on the first energy
landscape. For example if we want to start with 35% of the population
in minimum 3 and 65% in minimum 5 the I<STRING> looks as follows
"3=.35:5=.65" (floats must sum to 1).

=item B<-eq> I<INT>

If speficied, treekin simulation on the last landscape is performed
until the pecified stop time (recommended value 1000000).

=item B<-rs> I<STRING>

Suffix of the binary rates files produces by barriers (default: rates.bin)

=item B<-recover>

Recover from pre-computed (binary) eigenvectors and values (e.g. computed in
parallel) as computed by Treekin's C<--dumpE> option.  Speeds up repeated
simulations with different stop times etc.  For input file C<foo.rates.bin>,
eigenvectors and -values need to be stored in files C<foo.evecs.bin> and
C<foo.evals.bin>.

=back

=head1 SOURCE AVAILABITY

Source code for this distribution is available from the
L<ViennaRNA/BarMap GitHub
repository|https://github.com/ViennaRNA/BarMap>.

=head1 AUTHORS

Christoph Flamm, Michael T. Wolfinger

=head1 BUGS

Please send comments and bug reports to E<lt>mtw@tbi.univie.ac.at<gt>.

=cut

__END__
# my $a1 = [qw/a b c d e f g/];
# my $b1 = [qw/h i j k l m n/];
# my $c1 = [qw/o p q r s t u/];

# my $a2 = zip_with_natural_numbers($a1, natural_numbers());
# my $b2 = zip_with_natural_numbers($b1, natural_numbers());
# my $c2 = zip_with_natural_numbers($c1, natural_numbers());

# print Dumper($a2);
# print Dumper($b2);

# my $cols = to_cols($#$a2);
# $cols->($a2);
# print Dumper($cols->());
# $cols->($b2);
# $cols->($c2);
# print Dumper($cols->());

#my $a3 = unzip($a2);
#print Dumper($a3);
