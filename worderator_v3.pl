#!/usr/bin/perl
#    the above line is the SHEBANG
#
# worderator.pl input_pairs_file.tbl
#
#   combine all word pairs from a tabular file
#   to create a graph in order to reconstruct original sentences
#   from where the word pairs came from.
#
use strict;
use warnings;
#
# Dumping complex data structures
use Data::Dumper;
#
# Higher time resolution for sleep
use Time::HiRes;
#
# Capturing time points in the program execution
use Benchmark;
#
# Adding some command-line arguments processing
use Getopt::Long;

# INITIALIZE VARIABLES

my @exectime = ();

my %GRAPH = (
	'P2C'      => {},
	'C2P'      => {},
        'Nnodes'   => 0,
        'NODES'    => [],
        'TPARENTS' => [],
        'BCHILDS'  => []
	);
    
my ($verbose, $filename, $outputfile, $outfile, $dotformat, $kmer
    ) = (0, undef, undef, undef, 'png', 20);
#          # default: quiet and reuse input file name for output

my %VALID_DOTFORMATS = (
    'png' => 0,
    'jpg' => 0,
    'gif' => 0,
    'pdf' => 0 );

#
## PROCESSING COMMAND-LINE ARGUMENTS

GetOptions(
    "v|verbose"     => sub { $verbose = 1; }, # Setting levels of "verbosity"
    "d|debug"       => sub { $verbose = 2; }, # default $verbose is 0 that means "be quiet"
    "o|outfile=s"   => \$outfile,
    "t|dotformat=s" => \$dotformat,
    "k|kmer=i"      => \$kmer,
    "h|help"        => \&help                 # Calling help function, passed by reference
);

# processing filenames

$filename = shift @ARGV;
if (defined($outfile)) {
    $outputfile = $outfile ;
} else {
    ($outputfile = $filename) =~ s/\.[^\.]+$//;
};

# check the dotformat for valid images formats accepted by dot program

$dotformat = lc($dotformat);
exists($VALID_DOTFORMATS{$dotformat}) || do {
    print STDERR "### ERROR ### $dotformat is not a valid graphical format for dot program...\n",
                 "###           Using default PNG format... (see help information).\n";
    $dotformat = 'png'; 
};

# Add a check to ensure kmer is valid
if ($kmer <= 1) {
    die "## ERROR ## You must provide a valid k-mer size greater than 1 using -k (e.g., -k 3).\n";
}

#
# MAIN BLOCK

push @exectime, (new Benchmark);
print STDERR "#####\n##### RUNNING $0 PID[$$] ",
             "$ENV{'USER'} #####\n#####\n"
             if $verbose;

&read_from_input_file($verbose, $filename, \%GRAPH);

&search_top_to_bottom($verbose, $outputfile."_childs.tbl",
                      $GRAPH{'P2C'}, $GRAPH{'BCHILDS'});
&search_top_to_bottom($verbose, $outputfile."_parents.tbl",
                      $GRAPH{'C2P'}, $GRAPH{'TPARENTS'});

&traverse_graph_function($verbose, $outputfile."_sentences.tbl",
                         $GRAPH{'P2C'}, $GRAPH{'TPARENTS'});
    # recursion can be improved to count how many "strings" were produced
    # by passing a counter by reference, like $sentences here
  
push @exectime, (new Benchmark);

print STDERR "#####\n##### Time spent: ",
             timestr(timediff($exectime[ $#exectime ],
                              $exectime[($#exectime - 1) ])),
             " #####\n#####\n",
             "##### $0 HAS FINISHED #####\n#####\n"
             if $verbose;

exit(0);

# FUNCTIONS

sub help() {

    print STDERR <<'EOHelp';
#
# USAGE:
#
#   worderator.pl [options] input_pairs_file.tbl
#
# DESCRIPTION:
#
#   Combine all word pairs from a tabular file to create a graph
#   in order to reconstruct original sentences from where
#   the word pairs came from.
#   It produces several files, one storing the generated graph 
#   in GraphViz dot format, one storing the parent nodes
#   (topmost starting nodes), another with the child leaves
#   (bottom-most terminal nodes), and another one with
#   the longests sentences we can reconstruct following the paths
#   available on the graph.
#   The program also runs "dot" command to generate a PNG figure
#   displaying the whole graph obtained from input file.
#
# OPTIONS:
#
#   -o | --outfile "file.out"
#       Setting output filename prefix for all the outputs;
#       default is set to input filename without
#       the last extension.
#
#   -t | --dotformat "format"
#       Default dot graphical output is PNG,
#       you can choose one of the following formats:
#         png, jpg, gif, pdf
#
#   -v | --verbose
#       Reporting program execution (quiet by default).
#
#   -d | --debug
#       Extended report from internal data structures.
#       DO NOT USE IN PRODUCTION
#              (only with test sets while developing).
#
#   -h | --help
#       Print this help.
#
EOHelp

    exit(1);

} # help

sub read_from_input_file() {

    my ($verbose, $filename, $graph) = @_;
        # this takes first argument into the first variable
        # and the second argument into the second variable
    my %UNIQUE = ();

    print STDERR "# Reading FASTQ from $filename and generating k-mers (k=$kmer)...\n"
                 if $verbose;
    
    open(IFH, $filename) ||
        die("## ERROR ## Cannot open $filename...\n");

    my $reads_processed = 0;
    my $kmers_count = 0;

    # READ FASTQ (4 lines per record)
    while (my $header = <IFH>) {
        my $seq  = <IFH>;
        my $plus = <IFH>;
        my $qual = <IFH>;

        
        chomp($seq); # Remove newline from sequence

        # Validate sequence length vs kmer size
        my $len = length($seq);
        next if $len < $kmer;

        $reads_processed++;

        
        # SLIDING WINDOW to generate K-mers
        # We go from 0 up to Length - K
        for (my $i = 0; $i <= ($len - $kmer); $i++) {
            
            # Extract the full k-mer
            my $current_kmer = substr($seq, $i, $kmer);

            # In De Bruijn graphs:
            # Parent is the Prefix (length k-1)
            # Child is the Suffix (length k-1)
            my $node_len = $kmer;
            my $wordA = substr($current_kmer, 0, $node_len); # Prefix
            my $wordB = substr($current_kmer, 1, $node_len); # Suffix

            # --- GRAPH CONSTRUCTION (Same as original code) ---
            
            # PARENTS
            exists($graph->{'P2C'}{$wordA}) || ($graph->{'P2C'}{$wordA} = {});
            exists($graph->{'P2C'}{$wordA}{$wordB}) || ($graph->{'P2C'}{$wordA}{$wordB} = 0);
            $graph->{'P2C'}{$wordA}{$wordB}++;

            # CHILDS
            exists($graph->{'C2P'}{$wordB}) || ($graph->{'C2P'}{$wordB} = {});
            exists($graph->{'C2P'}{$wordB}{$wordA}) || ($graph->{'C2P'}{$wordB}{$wordA} = 0);
            $graph->{'C2P'}{$wordB}{$wordA}++;
            
            # Add to UNIQUE list
            $UNIQUE{$wordA}++;
            $UNIQUE{$wordB}++; 
            $kmers_count++;
        }; 

    }; # <IFH>
    close(IFH);

    print STDERR "## Processed $reads_processed reads, generated $kmers_count edges...\n" if $verbose;

    # Populate the node list
    foreach my $node (keys %UNIQUE) {
        push @{ $graph->{'NODES'} }, $node;
        $graph->{'Nnodes'}++;
    }

} # read_from_input_file



sub search_top_to_bottom() {
    my ($verbose, $outputfile, $graph, $array) = @_;

    print STDERR "# Writing node list to $outputfile...\n"
                 if $verbose;
    
    open(CHILDFH, "> $outputfile") ||
        die("## ERROR ## Cannot open $outputfile...\n");

    my $nc = 0;
    foreach my $parent (keys %{$graph}) {
        foreach my $child (keys %{$graph->{$parent}}) {
             (exists($graph->{$child})) || do {
                 push @$array, $child;
                 print CHILDFH $child,"\n";
                 $nc++;
             };
        }; # foreach $child
    }; # foreach $parents
    
    close(CHILDFH);
    
    print STDERR "# SAVED $nc terminal nodes...\n# NODES: @{$array}\n"
          if $verbose;

} # search_top_to_bottom

sub traverse_graph_function() {
    my ($verbose, $outputfile, $graph, $array) = @_;
    my $sentences = 0;
 
    print STDERR "# Retrieving sentences from graph to $outputfile...\n"
          if $verbose;
   
    open(SENTENCESFH, "> $outputfile") ||
        die("## ERROR ## Cannot open $outputfile...\n");

    # looping through parent topmost nodes
    foreach my $PARENT (@$array) {
    
        &recursion_on_the_graph($verbose, \*SENTENCESFH, $graph,
                                \$sentences, $PARENT, $PARENT);
    
    };
    
    close(SENTENCESFH);
    
    print STDERR "# Generation of ${sentences} sentences COMPLETED...\n"
                 if $verbose;

} # traverse_graph_function

sub recursion_on_the_graph() {

    # Disabling the warnings on "Deep recursion"
    no warnings 'recursion';

    # at each recursion level the @words array grows with another node
    # because we provide an extra argument to this function after this local array
    my ($verbose, $ofh, $graph, $sentences, $parent, @words) = @_;

    my $n = sprintf("%09d", scalar @words);
    
    foreach my $child (keys %{ $graph->{$parent} }) {
     
        # We have used the hash of hashes as a counter of
        # the times we have found a given pair of words
        # (see read_from_input_file function for the "++").
        $graph->{$parent}{$child}--;

        ## just to debug recursion show vars through STDERR
        print STDERR "$n : $parent -> $child : $graph->{$parent}{$child}\n"
              if $verbose==2;

        # first recursion exit condition:
        #     no more child (terminal node).
        # second recursion exit condition:
        #     the edge has been visited as many times as
        #     the number of times it was found from input data.
        # As we are counting the number of visits done on each edge,
        # we can consider this function an implementation of a simplified
        # Eulerian path search algorithm (yet not optimal of coursse).
        if (!exists($graph->{$child}) || $graph->{$parent}{$child} < 0) {
            print $ofh join(" ", @words, $child), "\n";
            $$sentences++; # we update string counter by dereferencing it
        } else {
            &recursion_on_the_graph($verbose, $ofh, $graph,
                                    $sentences, $child, @words, $child);
        };

    }; # foreach $child

} # recursion_on_the_graph

