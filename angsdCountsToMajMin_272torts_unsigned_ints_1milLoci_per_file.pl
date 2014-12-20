#!/usr/bin/perl


# For our unbiased estimates of pairwise divergence, we want to run calculations
# on hundreds of millions of sites, where we have read counts for the major and
# minor alleles.
#
# This script does a few things. It takes as input the read counts dumped by
# ANGSD, finds the minor and major read counts for each individual at that locus,
# then packs them into binary file representations of that data. The binary file
# format is as follows:
#    -Every read depth (integer) for the major or minor allele at
#     a locus is encoded as an 8-bit unsigned integer (so values can range from 0 to 255).
#     In this experiment, we don't want to include loci with read depth greater
#     than 255, because with coverage around 1X, such read depth would be an indication
#     that something is wrong. In fact, we will cap the maximum read depth to consider
#     a locus at a much lower number (10 or 15, set with $maxReadDepth)
#    -Every individual gets two integers in the output file for each locus, with
#     the major allele read count first, and the minor allele read count second
#        -As such, there are N*2 integers per locus, where N is the number of bams
#         (individuals) input into ANGSD, in the same order as the bam file list
#    -Since we eventually want to be able to process many loci in parallel, the
#     output will be split into several/many output files. Currently we are working
#     with a 32-CPU machine, so we will initially be splitting the output into 30
#     separate output files. If we assume the genome is approximately 1.9 gb,
#     that means we need a new output file every ~ 63 million loci.


use strict;
use warnings;
use List::Util qw(sum max);
use Data::Dumper;
use Getopt::Long;

my $MAF;
my $countsFile;
my $help;

GetOptions  ("maf=s"      => \$MAF,
             "counts=s"   => \$countsFile,
             "help|man"   => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$MAF or !$countsFile or $help) {
    die "Must supply --maf and --counts.\n";
}

my @etortNumbers = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,164,166,167,168,169,170,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,202,203,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,282,283,284,285,286,287,288,289,296,297);

my $maxReadDepth = 900;
my $locusCounter = 0;
my $fileCounter = 0;
my $fileName = "272torts_1e-6_minQ30." . $fileCounter . "mil.binary8bitunsigned";
my $outFH;

my %alleleOrderHash = ("A" => 0,
                       "C" => 1,
                       "G" => 2,
                       "T" => 3
                       );

open(my $MAFfh, "<", $MAF) or die "Couldn't open $MAF for reading: $!\n";
open(my $countsFH, "<", $countsFile) or die "Couldn't open $countsFile for reading: $!\n";
open($outFH, ">", $fileName) or die "Couldn't open $fileName for writing: $!\n";

while ( not eof $countsFH and not eof $MAFfh ) { # Read from STDIN, so we can "gunzip -c" right into it from the ANGSD output
    my $countsLine = <$countsFH>;
    my $MAFline = <$MAFfh>;
    chomp($countsLine);
    if ($countsLine =~ /^ind/) { # This is the header line of the counts file
        next;
    }

    if ($locusCounter % 1000000 == 0) {
        $fileCounter++;
        close($outFH);
        $fileName = "272torts_1e-6_minQ30." . $fileCounter . "mil.binary8bitunsigned"; 
        open($outFH, ">", $fileName) or die "Couldn't open $fileName for writing: $!\n";
    }

    my @mafFields = split(/\t/, $MAFline);
    my $majorAllele = $mafFields[2];
    my $minorAllele = $mafFields[3];
    
    my @countsFields = split(/\t/, $countsLine);
    if (sum(@countsFields) > $maxReadDepth) {
        next;
    }
    
    my $tortCounter = 0;
    while ($tortCounter < scalar(@etortNumbers)) {
        my @individualDepths = splice(@countsFields,0,4); # Remove the first 4 elements of the array    
        my $majorDepth = $individualDepths[$alleleOrderHash{$majorAllele}];
        my $minorDepth = $individualDepths[$alleleOrderHash{$minorAllele}];
        
        #print $outFH "$majorDepth\t$minorDepth\t";
        print $outFH pack("C", $majorDepth);  # "C" = 8-bit integer
        print $outFH pack("C", $minorDepth);   
        $tortCounter++;   
    }
    $locusCounter++;
}

