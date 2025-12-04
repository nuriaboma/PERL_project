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

my ($verbose, $filename, $outputfile, $outfile, $dotformat, $K) = (0, undef, undef, undef, 'png', 21);

my %VALID_DOTFORMATS = (
    'png' => 0,
    'jpg' => 0,
    'gif' => 0,
    'pdf' => 0 
);


# COMMAND-LINE OPTIONS
GetOptions(
    "v|verbose" => sub { $verbose = 1; },
    "d|debug"   => sub { $verbose = 2; },
    "o|outfile=s" => \$outfile,
    "t|dotformat=s" => \$dotformat,
    "k|kmer=i" => \$K,
    "h|help" => \&help
);

# processing filenames

$filename = shift @ARGV;
if (defined($outfile)) {
    $outputfile = $outfile ;
} else {
    ($outputfile = $filename) =~ s/\.[^\.]+$//;# if defined $filename;
};

$dotformat = lc($dotformat);
exists($VALID_DOTFORMATS{$dotformat}) || do {
    print STDERR "### ERROR ### $dotformat is not a valid graphical format for dot program...\n",
                 "###           Using default PNG format... (see help information).\n";
    $dotformat = 'png'; 
};

# Add a check to ensure kmer is valid
if ($K <= 1) {
    die "## ERROR ## You must provide a valid k-mer size greater than 1 using -k (e.g., -k 30).\n";
}


# MAIN
push @exectime, (new Benchmark);
print STDERR "##### RUNNING $0 PID[$$] #####\n" if $verbose;

read_from_input_file($verbose, $filename, \%GRAPH);

search_top_to_bottom($verbose, $outputfile."_childs.tbl", $GRAPH{'P2C'}, $GRAPH{'BCHILDS'});

search_top_to_bottom($verbose, $outputfile."_parents.tbl", $GRAPH{'C2P'}, $GRAPH{'TPARENTS'});

traverse_graph_function($verbose, $outputfile."_contigs.tbl", $GRAPH{'P2C'}, $GRAPH{'TPARENTS'});

print STDERR "#####\n##### Time spent: ",
             timestr(timediff($exectime[ $#exectime ],
                              $exectime[($#exectime - 1) ])),
             " #####\n#####\n",
             "##### $0 HAS FINISHED #####\n#####\n"
             if $verbose;

exit(0);

# ----------------- FUNCTIONS -----------------

sub help {
    print STDERR <<'EOHelp';
USAGE: 

    worderator.pl [options] input_reads.fastq

DESCRIPTION:

    Build a de-Bruijn graph from FASTQ reads using k-mers, list terminal nodes
    and parents, and output the longest contig sequence.

OPTIONS:

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
    die("## ERROR ## No input file provided!\n") unless defined $filename; #No need to try to open file if file was not provided

    print STDERR "# Reading FASTQ from $filename, building de-Bruijn graph (k=$K)...\n" if $verbose;

    open(IFH, $filename) ||
        die("## ERROR ## Cannot open $filename...\n");

    my $pairs=0; 
    my $n=0; 
    my $k=$K; 
    my $node_len=$k-1;

    while (1) {
        my $header = <IFH>; last unless defined $header;
        my $seq    = <IFH>;
        my $plus   = <IFH>;
        my $qual   = <IFH>; #last unless defined $qual;

        chomp($seq); 
        $n+=4;

        next if length($seq)<$k; #chack sequence length is less than k

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


sub traverse_graph_debruijn_dfs_collect {

    no warnings 'recursion';

    my ($verbose, $graph, $node, $sequence, $contigs_ref) = @_;

    foreach my $child (keys %{ $graph->{$node} }) {

        # Decrement edge usage counter
        $graph->{$node}{$child}--;

        print STDERR "DFS: $node -> $child (remaining: $graph->{$node}{$child})\n"
            if $verbose == 2;

        # Check terminal conditions
        if (!exists($graph->{$child}) || $graph->{$node}{$child} < 0) {

            my $final_seq = $sequence . substr($child, -1, 1);

            push @$contigs_ref, $final_seq;

            # print STDERR "  TERMINAL -> contig length=", length($final_seq), "\n"
            #     if $verbose;

        } else {

            my $next_seq = $sequence . substr($child, -1, 1);

            &traverse_graph_debruijn_dfs_collect(
                $verbose,
                $graph,
                $child,
                $next_seq,
                $contigs_ref
            );
        }
    }
}


sub traverse_graph_function {

    my ($verbose, $outputfile, $graph, $top_parents) = @_;

    print STDERR "# Traversing graph and collecting contigs...\n"
        if $verbose;

    my @all_contigs = ();

    # ---- DFS for each starting parent ----
    foreach my $parent (@$top_parents) {

        &traverse_graph_debruijn_dfs_collect(
            $verbose,
            $graph,
            $parent,
            $parent,      # initial sequence is the whole k-mer
            \@all_contigs
        );
    }

    print STDERR "# Total contigs generated: ", scalar @all_contigs, "\n"
        if $verbose;

    # ---- Sort contigs by length descending ----
    my @sorted = sort { length($b) <=> length($a) } @all_contigs;

    # ---- Keep only 5 longest ----
    my @top5 = @sorted[0 .. ($#sorted < 4 ? $#sorted : 4)];

    # ---- Output to file ----
    open(my $OUT, ">", $outputfile)
        or die "### ERROR: Cannot open $outputfile\n";

    my $count = 1;
    foreach my $seq (@top5) {
        my $len = length($seq);

        print $OUT "Sequence $count: $len bases\n";
        print $OUT "$seq\n\n";

        $count++;
    }

    close($OUT);

    print STDERR "# Wrote top 5 contigs to $outputfile\n";
}



# ---------------- memoized DFS -----------------
# my %LONGEST_PATH;

# sub get_longest_path_from_node {
#     no warnings 'recursion';
#     my ($graph, $node) = @_;
#     return $LONGEST_PATH{$node} if exists $LONGEST_PATH{$node};
#     return $node unless exists $graph->{$node};

#     my $best_seq='';
#     foreach my $child (keys %{$graph->{$node}}) {
#         my $child_seq = get_longest_path_from_node($graph, $child);
#         my $merged = $node . substr($child_seq,1);
#         $best_seq = $merged if length($merged)>length($best_seq);
#     }
#     return $LONGEST_PATH{$node} = $best_seq;
# }

# sub traverse_graph_function {

#     my ($verbose, $outputfile, $graph, $top_parents) = @_;
    
#     print STDERR "# Retrieving contigs from graph to $outputfile...\n"
#           if $verbose;

#     open(CONTIGSFH, ">", $outputfile) ||
#         die("## ERROR ## Cannot open $outputfile...\n");

#     my $longest='';
#     foreach my $parent (@$top_parents) {
#         my $seq = &get_longest_path_from_node($graph, $parent);
#         print CONTIGSFH $seq,"\n" if $verbose;
#         $longest = $seq if length($seq) > length($longest);
#     }
#     close(CONTIGSFH);
#     print STDERR "# Longest sequence length: ", length($longest), "\n";
# }
