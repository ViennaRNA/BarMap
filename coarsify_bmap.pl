#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2019-02-16 11:44:52 ivo>

# input at least barriers .bar file and a treekin trajectory
# output coarse grained dynamic obtained from merging lmins
#
# If given, regexes mark special features, only lmins with the same
# tags can be merged.  oreach lmin list which regex in @regs matched

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Path::Class;
use Carp;
use strict;
use warnings;

my $minh=3;
my $outdir = undef;
my $regsfile = undef;
my $bmfile = undef;
my $tkfile = undef;
my $debug_tags;
my @merge_lists = ();
my @new_labels = ();
my @regs = ();
my @tags = ("");

Getopt::Long::config("no_ignore_case");
pod2usage(-verbose => 1)
  unless GetOptions("minh|m=f"   => \$minh,
		    "outdir|o=s" => \$outdir,
		    "regs|r=s"   => \$regsfile,
		    "barmap|b=s" => \$bmfile,
		    "tkin|t=s"   => \$tkfile,
		    "debug-tags" => \$debug_tags,
		    "man"        => sub{pod2usage(-verbose => 2)},
		    "help|h"     => sub{pod2usage(-verbose => 1)}
		   )
  or pod2usage( "Try '$0 --help' for more information." );

unless (defined ($tkfile)){
  warn "Please provide an input (multi)treekin file via -t|--tkin option ...";
  pod2usage(-verbose => 1);
}

if(scalar(@ARGV)<1){
  warn "Please provide least one barriers file on the command line.";
  pod2usage(-verbose => 1 );
}

if($regsfile){
  open my $rfh, "<", $regsfile or
    croak "Could not open the specified regular expression file >",$regsfile,"<";
  foreach(<$rfh>){
    chomp;
    push(@regs, $_)
      if($_ !~ /^#/ and $_ !~ /^\s*$/); # ignore comment and empty lines
  }
  close $rfh;
}

unless (defined $outdir){
  $outdir = "coarse_$minh";
}
mkdir $outdir;
croak "cannot create directory $outdir" unless [-d $outdir];

# process bar file(s)
for my $nf (0..$#ARGV) {
  my $curfile = $ARGV[$nf];
  my ($seq, @lmin) = read_bar($curfile);
  @tags = ("");
  foreach (@lmin[1..$#lmin]) { # [0] is empty
    my $str = $_->[1];
    my $t = "";
    foreach my $r (@regs) {
      $t .= ($str =~ /$r/) ? '1' : '0';
    }
    push @tags, $t;
  }

  my @sbar = sort {$b->[4] <=> $a->[4]} @lmin[1..$#lmin];
  # print Dumper(\@sbar);
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

  my $bn = basename($ARGV[$nf], ".bar");
  my $newbar = file($outdir,$bn.".bar");
  #print "newbar ". Dumper("$newbar");
  open(my $barfh, ">", $newbar) or carp "can't open > $newbar";

  my @new = print_bar($barfh, $seq, \@lmin, \@merge);  # @new contains the new lmin labels
  close($barfh);
  push @merge_lists, \@merge;
  push @new_labels, \@new;
}

my $bntk = basename($tkfile);
my $newtkin = file($outdir,$bntk);
open(my $tkfh, ">", "$newtkin") or carp "can't open > $newtkin";
process_treekin($tkfh, $tkfile, \@merge_lists);
close($tkfh);

if($bmfile){
  my $bnbm = basename($bmfile);
  my $newbm = file($outdir,$bnbm);
  open(my $bmfh, ">", "$newbm") or carp "can't open > $newbm";
  process_barmap($bmfh, $bmfile, \@merge_lists, \@new_labels);
  close($bmfh);
}

sub read_bar {
  my ($barfile) = @_;
  open my $bf, "<", $barfile or carp "Cannot open barriers file for reading";
  $_ = <$bf>;
  warn "no seq in bar file $ARGV" unless /^\s+(\S+)/;
  my $seq = $1;
  my @lmin;
  my $nn = 1;
  while (<$bf>) {
    my @F = split;
    next if @F < 5; # shouldn't happen
    $lmin[$nn++] = \@F;
    last if eof;
  }
  close $bf;
  #print Dumper(\@lmin);
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
    printf $fh "%3d %s %6.2f %4d %6.2f", $n, $str, $e, $new[$f], $b;
    print $fh " $tags[$num]" if $debug_tags;
    print $fh "\n";

    $new[$num] = $n++;
  }
  return @new;
}

sub process_treekin {
  my ($fh, $tkfile, $merge_l) = @_;
  my $head;
  open my $tk, "<", $tkfile or carp "Cannot open barriers file for reading";
  foreach my $merge (@$merge_l) {
    while (<$tk>) {
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
    <$tk>;
  }
  close $tk;
}

sub process_barmap {
  my ($fh, $bf, $merge_lists, $new_labels) = @_;
 # print Dumper($merge_lists);
  # print Dumper($new_labels);
  open my $handle, "<", $bf or carp "Cannot open barmap file for reading";
  my @lines = <$handle>;
  close $handle;
  #print Dumper(\@lines);
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

coarsify_bmap.pl - coarse grain barriers, treekin and barmap files, for
better visualization and analysis

=head1 DESCRIPTION

Merge connected macro-states of a barriers files if the barrier height
is below the selected --minh value and/or states contain similar
structural elements, specified as regular expressions. If both
B<--minh|-m> and B<--regs|-s> is used macro-states that contain
different structural elements are never merged although the barrier
height might be lower than the specified minimum. The file containing
treekin trajectories (foo.mp) is re-written by summing up the
population densities of any merged states. A new barmap file,
representing the above merged macrostates is created if the
B<--barmap|-b> option is provided . All output files are written to
B<--outdir|-o> directory.

=head1 SYNOPSIS

  coarsify_bmap.pl [--minh minh] [--regs file] [--barmap foo.barmap]
                   [--tkin foo.mp] [--outdir out]
                   foo001.bar [foo002.bar foo003.bar [...]]

=head1 OPTIONS

=over 5

=item B<--minh|-m>

Local minima that are connected by a barrier hight lower than this
value are merged into one minimun/basin S<(default B<3>)>.

=item B<--regs|-r>

File that specifies structural elements as regular expressions
S<(default B<None>)>.

=item B<--barmap|-b>

Barmap output file to be processed S<(default B<None>)>.

=item B<--tkin|-t>

(Multi)treekin output file to be processed (B<required>). 

=item B<--outdir>

Specifiy the name of an output directory where all created files end
up S<(default B<coarse_minh>)>.

=item B<--debug-tags>

Append the matched tags to each line of the new bar file
in order to verify which minima match which regex.

=back

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

=head1 SOURCE AVAILABITY

Source code for this distribution is available from the
L<ViennaRNA/BarMap GitHub
repository|https://github.com/ViennaRNA/BarMap>.

=head1 AUTHORS

Ivo L Hofacker, Michael T Wolfinger, Sven Findeiss

=head1 BUGS

If in doubt our program is right, nature is at fault. Please send
comments and bug reports to E<lt>mtw@tbi.univie.ac.atE<gt>

=cut

__END__
