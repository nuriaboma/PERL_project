#!/usr/bin/perl
# Build a de-Bruijn graph from FASTQ reads and output longest contigs.

use strict;
use warnings;
use Data::Dumper;
use Time::HiRes;
use Benchmark;
use Getopt::Long;

# VARIABLES
my %GRAPH = (
    'P2C'      => {},   # parent -> child -> count
    'C2P'      => {},
    'TPARENTS' => [],   # top parents (start nodes)
    'BCHILDS'  => []    # bottom childs (terminal nodes)
);

my ($verbose, $filename, $outputfile, $outfile, $K) = (0, undef, undef, undef, 21);


# COMMAND-LINE OPTIONS
GetOptions(
    "v|verbose" => sub { $verbose = 1; },
    "d|debug"   => sub { $verbose = 2; },
    "o|outfile=s" => \$outfile,
    "k|kmer=i" => \$K,
    "h|help" => \&help
);

# processing filenames

$filename = shift @ARGV;
if (defined($outfile)) {
    $outputfile = $outfile ;
} else {
    ($outputfile = $filename) =~ s/\.[^\.]+$// if defined $filename;
};


# Add a check to ensure kmer is valid
if ($K <= 1) {
    die "## ERROR ## You must provide a valid k-mer size greater than 1 using -k (e.g., -k 30).\n";
}


# MAIN
my $t_start = Benchmark->new;
print STDERR "##### RUNNING $0 PID[$$] #####\n" if $verbose;

read_from_input_file($verbose, $filename, \%GRAPH);

search_top_to_bottom($verbose, $outputfile."_childs.tbl", $GRAPH{'P2C'}, $GRAPH{'BCHILDS'});

search_top_to_bottom($verbose, $outputfile."_parents.tbl", $GRAPH{'C2P'}, $GRAPH{'TPARENTS'});

traverse_graph_function($verbose, $outputfile."_contigs.tbl", $GRAPH{'P2C'}, $GRAPH{'TPARENTS'});

my $t_end = Benchmark->new;
print STDERR "\n## Time spent: ",
    timestr(timediff($t_end, $t_start)),
    "\n\n##### $0 HAS FINISHED #####\n"
    if $verbose;

exit(0);

# ----------------- FUNCTIONS -----------------

sub help {
    print STDERR <<'EOHelp';
USAGE: 

    de-bruijn.pl [options] input_reads.fastq

DESCRIPTION:

    Build a de-Bruijn graph from FASTQ reads using k-mers, list terminal nodes
    and parents, and output the longest contig sequence.

OPTIONS:

    -o | --outfile "file.out"      Output filename prefix
    -k | --kmer <int>              Set k-mer size (default 21)
    -v | --verbose                 Report execution
    -d | --debug                   Extended internal debug
    -h | --help                    Print this help

EOHelp
    exit(1);
}

sub read_from_input_file {
    my ($verbose, $filename, $graph) = @_;
    die("## ERROR ## No input file provided!\n") unless defined $filename; #No need to try to open file if file was not provided

    print STDERR "# Reading FASTQ from $filename, building de-Bruijn graph (k=$K)...\n" if $verbose;

    open(IFH, $filename) ||
        die("## ERROR ## Cannot open $filename...\n");

    my $pairs=0; 
    my $n=0; 
    my $k=$K; 
    my $node_len=$k-1;

    while (1) {
        my $header = <IFH>; last unless defined $header; # End loop if there is no header line
        my $seq    = <IFH>;
        my $plus   = <IFH>;
        my $qual   = <IFH>; 

        chomp($seq); 
        $n+=4;

        next if length($seq)<$k; #check sequence length is less than k

        for(my $i=0; $i<=length($seq)-$k; $i++){
            my $kmer = substr($seq,$i,$k);

            # In De Bruijn graphs:
            # Parent is the Prefix (length k-1)
            # Child is the Suffix (length k-1)
            my $prefix = substr($kmer, 0, $node_len);
            my $suffix = substr($kmer, 1, $node_len);

            $graph->{'P2C'}{$prefix}{$suffix}++;
            $graph->{'C2P'}{$suffix}{$prefix}++;

            $pairs++;
        }
    }

    close(IFH);
    print STDERR "## Processed $pairs k-mer edges from $n FASTQ lines\n" if $verbose;
}

sub search_top_to_bottom {
    my ($verbose, $outputfile, $graph, $array) = @_;

    print STDERR "# Writing node list to $outputfile...\n"
                 if $verbose;

    open(OFH, "> $outputfile") ||
        die("## ERROR ## Cannot open $outputfile...\n");

    my $nc=0;
    foreach my $parent (keys %{$graph}) {
        foreach my $child (keys %{$graph->{$parent}}) {
            (exists($graph->{$child})) || do {
                 push @$array, $child;
                 print OFH $child,"\n";
                 $nc++;
            };
        };
    };
    close(OFH);
    print STDERR "# SAVED $nc terminal nodes\n" if $verbose;
};

# ---------------------------
# Compact linear unitigs in the De Bruijn graph
# ---------------------------
sub compact_unitigs {
    my ($graph, $k) = @_;

    my (%in, %out);
    foreach my $u (keys %$graph) {
        foreach my $v (keys %{ $graph->{$u} }) {
            $out{$u}++;
            $in{$v}++;
        }
    }
    foreach my $u (keys %$graph) {
        $in{$u}  ||= 0;
        $out{$u} ||= 0;
    }

    my %visited;
    my %unitigs;      # UID => sequence
    my %new_graph;    # UID adjacency
    my %node2unitig;  # original node => UID

    my $unitig_id = 0;

    foreach my $node (keys %$graph) {
        next if $visited{$node};

        # terminal or branching node
        if ($in{$node} != 1 || $out{$node} != 1) {

            my $current = $node;
            my $seq     = $node;
            my @members = ($node);
            $visited{$node} = 1;

            # walk forward while linear
            while (($out{$current} // 0) == 1) {
                my ($next) = keys %{ $graph->{$current} // {} };  # safe if $current not in graph
                last if $visited{$next};
                last if ($in{$next} // 0) != 1;  # safe

                $seq .= substr($next, -1);
                push @members, $next;

                $visited{$next} = 1;
                $current = $next;
            }

            # store unitig
            my $uid = "U$unitig_id";
            $unitig_id++;
            $unitigs{$uid} = $seq;
            foreach my $n (@members) {
                $node2unitig{$n} = $uid;
            }

            # add outgoing edges
            if (exists $graph->{$current}) {
                foreach my $child (keys %{ $graph->{$current} }) {
                    my $target_uid = $node2unitig{$child} // "U$unitig_id";
                    $new_graph{$uid}{$target_uid}++;
                }
            }
        }
    }

    return (\%new_graph, \%unitigs, \%node2unitig);
}

# ---------------------------
# Map node to unitig UID
# ---------------------------
sub find_unitig_for_node {
    my ($node, $node2unitig) = @_;
    return exists $node2unitig->{$node} ? $node2unitig->{$node} : undef;
}

# ---------------------------
# Eulerian traversal
# ---------------------------
sub eulerian_path_from_node {
    my ($start, $graph) = @_;
    my @stack = ($start);
    my @path;

    # local copy to avoid destroying original
    my %local = map { $_ => { %{$graph->{$_}} } } keys %$graph;

    while (@stack) {
        my $v = $stack[-1];

        if ($local{$v} && keys %{$local{$v}}) {
            my ($u) = keys %{$local{$v}};
            $local{$v}{$u}--;
            delete $local{$v}{$u} if $local{$v}{$u} == 0;
            push @stack, $u;
        } else {
            push @path, pop @stack;
        }
    }
    return @path;   # reversed
}

# ---------------------------
# Convert Euler path (unitigs) to full sequence
# ---------------------------
sub unitig_path_to_sequence {
    my ($path_ref, $unitigs) = @_;
    my @path = reverse @$path_ref;

    return undef unless @path;

    my $seq = $unitigs->{$path[0]};
    for my $i (1 .. $#path) {
        my $u       = $path[$i];
        my $unitig  = $unitigs->{$u};
        next unless defined $unitig;

        # append full unitig except first base
        my $append = substr($unitig, 1);
        $seq .= $append;
    }
    return $seq;
}

# ---------------------------
# Main traversal and top-5 contigs output
# ---------------------------
sub traverse_graph_function {
    my ($verbose, $outputfile, $graph, $top_parents, $k) = @_;

    print STDERR "# Compacting unitigs...\n" if $verbose;
    my ($compact_graph, $unitigs, $node2unitig) = &compact_unitigs($graph, $k);

    print STDERR "# Unitigs: ", scalar(keys %$unitigs), "\n" if $verbose;
    print STDERR "# Euler traversal on compacted graph...\n" if $verbose;

    my @contigs;

    foreach my $start (@$top_parents) {
        my $start_uid = &find_unitig_for_node($start, $node2unitig);
        next unless defined $start_uid;

        my @euler = &eulerian_path_from_node($start_uid, $compact_graph);
        my $seq   = &unitig_path_to_sequence(\@euler, $unitigs);

        push @contigs, $seq if defined $seq && length($seq) > 0;
    }

    # sort descending by length
    @contigs = sort { length($b) <=> length($a) } @contigs;

    my @top5 = @contigs[0 .. ($#contigs < 4 ? $#contigs : 4)];

    open(my $OUT, ">", $outputfile) or die "$!";
    my $n = 1;
    foreach my $seq (@top5) {
        print $OUT "Sequence $n: ", length($seq), " bases\n$seq\n\n";
        $n++;
    }
    close($OUT);

    print STDERR "# Finished. Top ", scalar(@top5), " sequences written to $outputfile\n" if $verbose;
};
