#!/usr/bin/env perl
# -*-Perl-*-
# Last changed Time-stamp: <2018-01-04 14:39:42 mtw>

# input barriers .bar file and a treekin trajectory
# output coarse grained dynamic obtained from merging lmins

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use strict;
use warnings;

my $minh=3;
my $outdir;

Getopt::Long::config("no_ignore_case");
pod2usage(-verbose => 1)
  unless GetOptions("minh=f"   => \$minh,
		    "outdir=s" => \$outdir,
		    "man"      => sub{pod2usage(-verbobose => 2)},
		    "help"     => sub{pod2usage(-verbose => 1)});

my @regs = (
	   );     # list of regexes that mark special features, only
                  # lmins with the same tags can be merged.
my @tags = ("");  # foreach lmin list which regex in @regs matched

$outdir = "coarse_$minh" unless defined($outdir);

mkdir $outdir;
die "can't creat directory $outdir"  unless [-d $outdir];

my @merge_lists;
my @new_labels;
# process the bar files
for my $nf (0..$#ARGV-2) {
  my ($seq, @lmin) = read_bar();

  foreach (@lmin) {
    my $str = $_->[1];
    my $t = "";
    foreach my $r (@regs) {
      $t .= ($str =~ /$r/) ? '1' : '0';
    }
    push @tags, $t;
  }

  my @sbar = sort {$b->[4] <=> $a->[4]} @lmin[1..$#lmin];
#  print Dumper(\@sbar);
  my @merge = ();
  $merge[$_] =$_ for (0..$#lmin);

  foreach my $s (@sbar) {
    next if $s->[4]>=$minh;
    next if $s->[0]==1; # just in case: never merge the ground state;
    my $f = $s->[3];
    next unless $tags[$s->[0]] eq $tags[$f];
    $f = $merge[$f] while ($merge[$f] != $f);
    $merge[$s->[0]] = $f;
  }

  my $newbar = "$outdir/$ARGV";
  open(my $fh, ">", $newbar) or die "can't open > $newbar";

  my @new = print_bar($fh, $seq, \@lmin, \@merge);  # @new contains the new lmin labels
  close($fh);
  push @merge_lists, \@merge;
  push @new_labels, \@new;
}
#print Dumper(\@merge_lists);

#foreach my $m (1..14) {
#  foreach my $l  (0..$#merge_lists) {
#    my $n = $merge_lists[$l][$m];
#    printf "%3d%2d ", $n, $new_labels[$l][$n] if defined($n);
#    print "      " if !defined($n);
#  }
#  print "\n";
#}


my $newtkin = $ARGV[-2];
open(my $fh, ">", "$outdir/$newtkin") or die "can't open > $outdir/$newtkin";
process_treekin($fh, \@merge_lists);
close($fh);

my $newbm = $ARGV[-1];
open(my $fhbm, ">", "$outdir/$newbm") or die "can't open > $outdir/$newbm";
process_barmap($fhbm, \@merge_lists, \@new_labels);
close($fhbm);


#---
sub read_bar {
  $_ = <>;
  warn "no seq in bar file $ARGV" unless /^\s+(\S+)/;
  my $seq = $1;
  my @lmin;
  my $nn = 1;
  while (<>) {
    my @F = split;
    next if @F < 5; # shouldn't happen
    #    splice(@F,2,1) if ($F[2] eq '('); # got 2 fields from e.g. "( -1.00)"
    #    $F[2] =~ s/[()]//g;               # remove brackets from energy field
    $lmin[$nn++] = \@F;
    last if eof;
  }
  return $seq, @lmin;
}

sub print_bar {
  my ($fh, $seq, $lmin, $merge) = @_;
  print $fh "    $seq\n";
  my $n=1;
  my @new = (0);
  foreach my $l (@{$lmin}[1..$#$lmin]) {
    my ($num, $str, $e, $f, $b) = @{$l};
    #    next if $f==0; # skip unconnected
    next unless $num == $merge->[$num];
    $f = $merge->[$f] while ($merge->[$f] < $f);
    printf $fh "%3d %s %6.2f $f%4d %6.2f\n", $n, $str, $e, $new[$f], $b;
    $new[$num] = $n++;
  }
  return @new;
}

sub process_treekin {
  my ($fh, $merge_l) = @_;
  my $head;
  foreach my $merge (@$merge_l) {
    while (<>) {
      print($fh $_) , last if /^&/;
      next if /^@/;    # get rid of xmgrace lines 
      print($fh $_) , next if /^#/; 
      print $fh "# coarse grained with min height $minh and regexes ",
	join("\t", @regs),"\n" unless $head;
      $head=1;   # done reading headlines, a "&" will mark the end of this trajectory
      my @l = split;
      my @ll = @l;
      
      printf $fh "%22.20e ", $l[0]; # print time
      for (my $i=$#l; $i>0; $i--) { # from right to left since we always have merge[l]<l;
	next unless $merge->[$i];
	next if $merge->[$i] == $i;
	$l[$merge->[$i]] += $l[$i];
      }
      my ($s1,$s2) = (0,0);
      for my $i (1..$#l) {
	$s2 += $ll[$i];
	next if (defined $merge->[$i] and $merge->[$i] != $i); # if I use absorbance I get more states than minima
	next if $merge->[$i] != $i;
	printf $fh "%e ", $l[$i];
	$s1 += $l[$i];
      }
      print $fh "\n";
      # print STDERR $s1-1, " ", $s2-1, " ", 1-$s1/$s2, "\n";
    }
    undef $head;
  }
  while (!eof) { # discard rest of file
    <>;
  }
}

sub process_barmap {
  my ($fh, $merge_lists, $new_labels) = @_;
  
  my @lines = <>;  # slurp barmap file
  print $fh shift(@lines);
  my @bmap;
  
  foreach (@lines) {
    my @row = map {s/ //g; $_ eq "" ? -1 : $_ } (m/(....)(?:...)?/g);
    push @bmap, [@row];
  }
  for my $l (0..$#$merge_lists) {
    my @merge = @{$merge_lists->[$l]};
    for my $r (0..$#bmap) {
      $bmap[$r][$l] = $new_labels->[$l][$merge[$bmap[$r][$l]]] 	if ($bmap[$r][$l]>0);
    }
  }
  foreach my $r (@bmap) {
    if ($r->[0]>0) {
      printf $fh "%3d", $r->[0];
    } else {
      print $fh "   ";
    }
    foreach my $m (@$r[1..$#$r]) {
      if ($m>0) {
	printf $fh " ->%3d", $m;
      } else {
	print $fh "      ";
      }
    }
    print $fh "\n";
  }
}

__END__

=head1 NAME

coarsify_bmap.pl - coarse grain barmap trajectories for easier
visualization and analysis

=head1 SYNOPSIS

coarsify_bmap.pl [--outdir out] [--minh minh] foo001.bar foo002.bar
[foo003.bar [...]] foo.mp foo.barmap

=head1 DESCRIPTION

coarse grain each bar file using the specified minimum barrier height
and/or regex. Re-write the file containing the treekin trajectories
(foo.mp) by summing up any merged states and create a new barmap file.
The new bar files, treekin file, and barmap file are written to the
outdir directory.

=head1 PUBLICATIONS

If you find this tool useful for your work or if you use any of its
components, please cite the following publications:

 BarMap: RNA Folding on Dynamic Energy Landscapes
 Ivo L. Hofacker, Christoph Flamm, Christian Heine, Michael T. Wolfinger,
 Gerik Scheuermann and Peter F. Stadler
 RNA (2010) 16:1308-1316 doi:10.1261/rna.2093310

 Efficient computation of RNA-ligand interaction dynamics
 Michael T. Wolfinger, Christoph Flamm, Ivo L. Hofacker
 Methods (2018)

=head1 AUTHORS

Ivo L Hofacker, Michael T Wolfinger

=head1 BUGS

If in doubt our program is right, nature is at fault. Please send
comments and bug reports to E<lt>mtw@tbi.univie.ac.atE<gt>

=cut

__END__
