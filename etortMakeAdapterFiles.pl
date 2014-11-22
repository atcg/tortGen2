#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $help = 0;
my $inFile;
my $outDir;

GetOptions  ("in=s"      => \$inFile,
             "out=s"      => \$outDir,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outDir or $help) {
    die "Must supply --in and --out.\n";
}

unless(-d $outDir) {
    mkdir $outDir;
}


open(my $inFH, "<", $inFile) or die "Couldn't open $inFile for reading: $!\n";

while (my $line = <$inFH>) {
    my @fields = split(/\t/, $line); # $fields[0] is the index sequence, $fields[1] is the sample name
    if ($fields[1] =~ /(.*)_R1.fastq.gz/) {
        my $outFile = $outDir . "/$1" . "_adapters.fasta";
        my $indexedAdapterSeq = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" . $fields[0] . "ATCTCGTATGCCGTCTTCTGCTTG";
        my $indexedAdapterRevComp = revcomp($indexedAdapterSeq);
        my $universalAdapterSeq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
        my $universalAdapterSeqRevComp = revcomp($universalAdapterSeq);
        
        open(my $outFH, ">", ">$outFile") or die "Couldn't open $outFile for writing: $!\n";
        print $outFH ">$1\_indexed_adapter\n$indexedAdapterSeq\n";
        print $outFH ">$1\_indexed_adapter_revcomp\n$indexedAdapterRevComp\n";
        print $outFH ">$1\_universal_adapter\n$universalAdapterSeq\n";
        print $outFH ">$1\_universal_adapter_revcomp\n$universalAdapterSeqRevComp\n";
    } else {
        die "Line should end with _R1.fastq.gz.\n";
    }
}







sub revcomp {
    my $dna = shift;

	# reverse the DNA sequence
    my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}