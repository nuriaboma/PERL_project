#!/usr/bin/perl
# worderator.pl input_pairs_file.fastq
# Build a de-Bruijn graph from FASTQ reads and output longest contigs.

use strict;
use warnings;
use Data::Dumper;
use Time::HiRes;
use Benchmark;
use Getopt::Long;

# VARIABLES
my @exectime = ();

my %GRAPH = (
    'P2C'      => {},   # parent -> child -> count
    'C2P'      => {},
    'TPARENTS' => [],   # top parents (start nodes)
    'BCHILDS'  => []    # bottom childs (terminal nodes)
);

my ($verbose, $filename, $outputfile, $outfile, $dotformat) = (0, undef, undef, undef, 'png');
my $K = 21;  # default k-mer size

my %VALID_DOTFORMATS = ( 'png'=>0,'jpg'=>0,'gif'=>0,'pdf'=>0 );

# COMMAND-LINE OPTIONS
GetOptions(
    "v|verbose" => sub { $verbose = 1; },
    "d|debug"   => sub { $verbose = 2; },
    "o|outfile=s" => \$outfile,
    "t|dotformat=s" => \$dotformat,
    "k|kmer=i" => \$K,
    "h|help" => \&help
);

$filename = shift @ARGV;
if (defined($outfile)) {
    $outputfile = $outfile ;
} else {
    ($outputfile = $filename) =~ s/\.[^\.]+$// if defined $filename;
}

$dotformat = lc($dotformat);
exists($VALID_DOTFORMATS{$dotformat}) or $dotformat = 'png';

# MAIN
push @exectime, (new Benchmark);
print STDERR "##### RUNNING $0 PID[$$] #####\n" if $verbose;

read_from_input_file($verbose, $filename, \%GRAPH);

# optional: skip for very large graphs
# save_dot_file($verbose, $outputfile, $GRAPH{'P2C'});
# run_graphviz_dot($verbose, $dotformat, $outputfile);

search_top_to_bottom($verbose, $outputfile."_childs.tbl", $GRAPH{'P2C'}, $GRAPH{'BCHILDS'});
search_top_to_bottom($verbose, $outputfile."_parents.tbl", $GRAPH{'C2P'}, $GRAPH{'TPARENTS'});

traverse_graph_function($verbose, $outputfile."_sentences.tbl", $GRAPH{'P2C'}, $GRAPH{'TPARENTS'});

push @exectime, (new Benchmark);
print STDERR "##### Time spent: ", timestr(timediff($exectime[-1], $exectime[-2])), " #####\n";
print STDERR "##### $0 HAS FINISHED #####\n";

exit(0);

# ----------------- FUNCTIONS -----------------

sub help {
    print STDERR <<'EOHelp';
USAGE: worderator.pl [options] input_reads.fastq

Build a de-Bruijn graph from FASTQ reads using k-mers, list terminal nodes
and parents, and output the longest contig sequence.

Options:
  -o | --outfile "file.out"      Output filename prefix
  -t | --dotformat "format"      PNG, JPG, GIF, PDF (default PNG)
  -k | --kmer <int>              Set k-mer size (default 21)
  -v | --verbose                 Report execution
  -d | --debug                   Extended internal debug
  -h | --help                    Print this help
EOHelp
    exit(1);
}

sub read_from_input_file {
    my ($verbose, $filename, $graph) = @_;
    die("## ERROR ## No input file provided!\n") unless defined $filename;

    print STDERR "# Reading FASTQ from $filename, building de-Bruijn graph (k=$K)...\n" if $verbose;

    open(my $IFH, $filename) or die "Cannot open $filename\n";

    my $pairs=0; my $n=0; my $k=$K; my $node_len=$k-1;
    die("## ERROR ## k must be >=2\n") if $node_len<1;

    while (1) {
        my $header = <$IFH>; last unless defined $header;
        my $seq    = <$IFH>;
        my $plus   = <$IFH>;
        my $qual   = <$IFH>; last unless defined $qual;

        chomp($seq); $n+=4;

        next if length($seq)<$k;

        for(my $i=0; $i<=length($seq)-$k; $i++){
            my $kmer = substr($seq,$i,$k);
            my $prefix = substr($kmer,0,$node_len);
            my $suffix = substr($kmer,1,$node_len);

            $graph->{'P2C'}{$prefix}{$suffix}++;
            $graph->{'C2P'}{$suffix}{$prefix}++;

            $pairs++;
        }
    }

    close($IFH);
    print STDERR "## Processed $pairs k-mer edges from $n FASTQ lines\n" if $verbose;
}

sub search_top_to_bottom {
    my ($verbose, $outputfile, $graph, $array) = @_;
    open(my $OFH, ">", $outputfile) or die "Cannot open $outputfile\n";

    my $nc=0;
    foreach my $parent (keys %{$graph}) {
        foreach my $child (keys %{$graph->{$parent}}) {
            push @$array, $child unless exists $graph->{$child};
            print $OFH $child,"\n" unless exists $graph->{$child};
            $nc++;
        }
    }
    close($OFH);
    print STDERR "# SAVED $nc terminal nodes\n" if $verbose;
}

# ---------------- memoized DFS -----------------
my %LONGEST_PATH;

sub get_longest_path_from_node {
    no warnings 'recursion';
    my ($graph, $node) = @_;
    return $LONGEST_PATH{$node} if exists $LONGEST_PATH{$node};
    return $node unless exists $graph->{$node};

    my $best_seq='';
    foreach my $child (keys %{$graph->{$node}}) {
        my $child_seq = get_longest_path_from_node($graph, $child);
        my $merged = $node . substr($child_seq,1);
        $best_seq = $merged if length($merged)>length($best_seq);
    }
    return $LONGEST_PATH{$node} = $best_seq;
}

sub traverse_graph_function {
    no warnings 'recursion';
    my ($verbose, $outputfile, $graph, $top_parents) = @_;

    open(my $OFH, ">", $outputfile) or die "Cannot open $outputfile\n";
    my $longest='';
    foreach my $parent (@$top_parents) {
        my $seq = get_longest_path_from_node($graph, $parent);
        print $OFH $seq,"\n" if $verbose;
        $longest = $seq if length($seq) > length($longest);
    }
    close($OFH);
    print STDERR "# Longest sequence length: ", length($longest), "\n";
}
