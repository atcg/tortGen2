#!/usr/bin/perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Cwd;
use Parallel::ForkManager;

my $help = 0;
my $readsDir;
my $adaptersDir;
my $outDir;
my $logFile;
my $threadsMax = 4;
my $trim;
my $filter; # If set, it will remove reads that were filtered by CASAVA
my $join;


GetOptions  ("reads=s"         => \$readsDir,
             "adapters=s"      => \$adaptersDir,
             "out=s"           => \$outDir,
             "log=s"           => \$logFile,
             "threads=i"       => \$threadsMax,
             "trim"            => \$trim,
             "filter"          => \$filter,
             "help|man" => \$help) || pod2usage(2);

if (!$readsDir or !$adaptersDir or !$outDir or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

my @samples = ("etort-1","etort-2","etort-3","etort-4","etort-5","etort-6","etort-7","etort-8","etort-9","etort-10","etort-11","etort-12","etort-13","etort-14","etort-15","etort-16","etort-17","etort-18","etort-19","etort-20","etort-21","etort-22","etort-23","etort-24","etort-25","etort-26","etort-27","etort-28","etort-30","etort-31","etort-32","etort-33","etort-34","etort-35","etort-36","etort-37","etort-38","etort-39","etort-40","etort-41","etort-42","etort-43","etort-44","etort-45","etort-46","etort-47","etort-48","etort-49","etort-50","etort-51","etort-52","etort-53","etort-54","etort-55","etort-56","etort-57","etort-58","etort-59","etort-60","etort-61","etort-62","etort-63","etort-64","etort-65","etort-66","etort-67","etort-68","etort-69","etort-70","etort-71","etort-72","etort-73","etort-74","etort-75","etort-76","etort-77","etort-78","etort-79","etort-80","etort-81","etort-82","etort-83","etort-84","etort-85","etort-86","etort-87","etort-88","etort-89","etort-90","etort-91","etort-92","etort-93","etort-94","etort-95","etort-96","etort-97","etort-98","etort-99","etort-100","etort-102","etort-103","etort-104","etort-105","etort-106","etort-107","etort-108","etort-109","etort-110","etort-111","etort-112","etort-113","etort-115","etort-116","etort-117","etort-119","etort-120","etort-121","etort-122","etort-123","etort-124","etort-125","etort-126","etort-127","etort-128","etort-129","etort-130","etort-131","etort-132","etort-133","etort-134","etort-135","etort-136","etort-137","etort-138","etort-139","etort-140","etort-141","etort-142","etort-143","etort-144","etort-145","etort-146","etort-147","etort-148","etort-149","etort-150","etort-151","etort-152","etort-153","etort-154","etort-155","etort-156","etort-157","etort-158","etort-159","etort-161","etort-162","etort-163","etort-164","etort-166","etort-167","etort-168","etort-169","etort-170","etort-172","etort-173","etort-174","etort-175","etort-176","etort-177","etort-178","etort-179","etort-180","etort-181","etort-182","etort-183","etort-184","etort-185","etort-186","etort-187","etort-188","etort-189","etort-190","etort-191","etort-192","etort-193","etort-194","etort-195","etort-196","etort-197","etort-202","etort-203","etort-208","etort-209","etort-210","etort-211","etort-212","etort-213","etort-214","etort-215","etort-216","etort-217","etort-218","etort-219","etort-220","etort-221","etort-222","etort-223","etort-224","etort-225","etort-226","etort-227","etort-228","etort-229","etort-230","etort-231","etort-232","etort-233","etort-234","etort-235","etort-236","etort-237","etort-238","etort-239","etort-240","etort-241","etort-242","etort-243","etort-244","etort-245","etort-247","etort-248","etort-249","etort-250","etort-251","etort-252","etort-253","etort-254","etort-255","etort-256","etort-257","etort-258","etort-259","etort-260","etort-261","etort-262","etort-263","etort-264","etort-265","etort-266","etort-267","etort-268","etort-269","etort-270","etort-271","etort-272","etort-273","etort-274","etort-275","etort-276","etort-277","etort-278","etort-282","etort-283","etort-284","etort-285","etort-286","etort-287","etort-288","etort-289","etort-296","etort-297","etort-298","etort-299","etort-300","etort-301");
open(my $logFH, ">", $logFile) or die "Couldn't open log file $logFile for writing: $!\n";


# Make sure we have the right directories for output files
my $startingDir = getcwd();
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make output directory $outDir: $!\n";
}

my $fqjDir = $outDir . "/fastq-join";
unless (-d $fqjDir) {
    mkdir $fqjDir;
}

my $trimmomaticDir = $outDir . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}



########## Remove the :Y: files that were designated as bad by CASAVA ##########
if ($filter) { 
    print $logFH "Removing reads that were filtered by CASAVA\n";
    foreach my $tort (@samples) {
        my $R1File = $readsDir . $tort . "_R1.fastq.gz";
        my $R2File = $readsDir . $tort . "_R2.fastq.gz";
        my $R1outFile = $readsDir . $tort . "R1_Ns.fastq.gz";
        my $R2outFile = $readsDir . $tort . "R2_Ns.fastq.gz";
        system("gunzip -c $R1File | fastq_illumina_filter -vN | gzip > $R1outFile");
        system("gunzip -c $R2File | fastq_illumina_filter -vN | gzip > $R2outFile");
    }
    print $logFH "Finished removing reads that were filtered by CASAVA\n\n";
}


########## Generate Trimmomatic commands and run Trimmomatic ##########

print $logFH "Generating Trimmomatic commands for all read files\n";
my @trimmomaticCommands;
foreach my $tort (@samples) {
    my $R1File = $readsDir . "$tort" . "_R1_Ns.fastq.gz";
    my $R2File = $readsDir . "$tort" . "_R2_Ns.fastq.gz";
    my $adaptersFile = $adaptersDir . $tort . "_adapters.fasta";
    my $R1OutFilePaired = $trimmomaticDir . "/$tort" . "_R1p_Ns_trim.fastq.gz";
    my $R2OutFilePaired = $trimmomaticDir . "/$tort" . "_R2p_Ns_trim.fastq.gz";
    my $R1OutFileSingles = $trimmomaticDir . "/$tort" . "_R1s_Ns_trim.fastq.gz";
    my $R2OutFileSingles = $trimmomaticDir . "/$tort" . "_R2s_Ns_trim.fastq.gz";

    push (@trimmomaticCommands, "java -Xmx8g -jar /home/evan/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 2 -phred33 $R1File $R2File $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adaptersFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");    
}
print $logFH "Finished generating trimmomatic commands for all tortoises\n";



# Now run Trimmomatic
if ($trim) {
    print $logFH "Running all trimmomatic commands\n";
    my $counter = 0;
    my $trimThreads = $threadsMax / 2; # We're using two threads for every Trimmomatic process
    my $trimForkManager = new Parallel::ForkManager($threadsMax);
    foreach my $trimCommand (@trimmomaticCommands) {
        $counter++;
        print $logFH "--------------------------------------------------\n";
        print $logFH "Trimmomatic command $counter: \n\t"; # Indent the next line to make it easier to find the commands in the text
        print $logFH $trimCommand . "\n";
        # sleep 10; # Turn this off because we're not actually running the commands so we don't need to be careful about jumbling output
        print "\n";
        $trimForkManager->start and next;
        print "\n";
        # system("$trimCommand"); # turn this off for now--starting the process at fastq-join
        print "Finished running the following:\n\t$trimCommand\n\n";
        $trimForkManager->finish;
    }
    $trimForkManager->wait_all_children;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running all trimmomatic commands\n";
    print $logFH "--------------------------------------------------\n\n";
}


# Now run fastq-join to merge overlapping reads
if ($join) {
    

}


    print $logFH "Running fastq-join to merge overlapping paired end reads\n";
    foreach my $readGroup (sort keys %sampleNamesHash) {
        sleep 10;
        print $logFH "\n\nStarting to process $readGroup through the fastq-join and assembly stages\n\n";
        my $combinedTrimmomaticSingles = $fqjDir . "/" . $readGroup . "_R1andR2trimmomaticSingles.fastq";
        my $R1pairedTrimmedUnzipped;
        my $R2pairedTrimmedUnzipped;
        if ($sampleNamesHash{$readGroup}{'R1_paired_trimmed'} =~ /(.*).gz$/) {
            $R1pairedTrimmedUnzipped = $1;
            system("gunzip $sampleNamesHash{$readGroup}{'R1_paired_trimmed'}");
        } else {
            $R1pairedTrimmedUnzipped = $sampleNamesHash{$readGroup}{'R1_paired_trimmed'}; # If it doesn't end in gz, it's already unzipped
        }
        if ($sampleNamesHash{$readGroup}{'R2_paired_trimmed'} =~ /(.*).gz$/) {
            $R2pairedTrimmedUnzipped = $1;
            system("gunzip $sampleNamesHash{$readGroup}{'R2_paired_trimmed'}");
        } else {
            $R2pairedTrimmedUnzipped = $sampleNamesHash{$readGroup}{'R2_paired_trimmed'}; # If it doesn't end in gz, it's already unzipped
        }
        my $fqjOutputPrefix = $fqjDir . "/" . $readGroup . "_trimmed_fqj.%.fastq";
        system("fastq-join -v ' ' $R1pairedTrimmedUnzipped $R2pairedTrimmedUnzipped -o $fqjOutputPrefix");
        system("cat $sampleNamesHash{$readGroup}{'R1_singles_trimmed'} $sampleNamesHash{$readGroup}{'R2_singles_trimmed'} > $combinedTrimmomaticSingles");
        #system("gunzip $combinedTrimmomaticSingles");
        my $joinedReads = $fqjDir . "/" . $readGroup . "_trimmed_fqj.join.fastq";
        my $joinedAndSingles = $fqjDir . "/" . $readGroup . "_combinedJoinedAndSingles.fastq";

        # Need to update this file name since we unzipped the file above 
        if ($combinedTrimmomaticSingles =~ /(.*).gz$/) {
            $combinedTrimmomaticSingles = $1;
        }
        
        system("cat $joinedReads $combinedTrimmomaticSingles > $joinedAndSingles");
    }
    
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running fastq-join on all samples\n";
    print $logFH "--------------------------------------------------\n\n";

}

if ($map) {
    my $mappingDir = $outDir . "/mapping";
    unless (-d $mappingDir) {
        mkdir $mappingDir;
    }
    foreach my $readGroup (sort keys %sampleNamesHash) {
        #my $singlesSamFile = $mappingDir . "/" . $readGroup . ".singlesAndJoined.sam";
        my $singlesBamFile = $mappingDir . "/" . $readGroup . ".singlesAndJoined.bam";
        #my $pairedSamFile = $mappingDir . "/" . $readGroup . ".paired.sam";
        my $pairedBamFile = $mappingDir . "/" . $readGroup . ".paired.bam";
        my $mergedBamFile = $mappingDir . "/" . $readGroup . "_merged.bam";
        my $reads1 = $fqjDir . "/" . $readGroup . "_trimmed_fqj.un1.fastq";
        my $reads2 = $fqjDir . "/" . $readGroup . "_trimmed_fqj.un2.fastq";
        my $readsSingles = $fqjDir . "/" . $readGroup . "_combinedJoinedAndSingles.fastq";
        system("bwa mem -t $threadsMax -M $reference $readsSingles | samtools view -@ $threadsMax -bS - > $singlesBamFile");
        system("bwa mem -t $threadsMax -M $reference $reads1 $reads2 | samtools view -@ $threadsMax -bS - > $pairedBamFile");
        
        #system("samtools view -@ 8 -bS $singlesSamFile > $singlesBamFile");
        #system("samtools view -@ 8 -bS $pairedSamFile > $pairedBamFile");
        
        system("samtools merge -@ $threadsMax $mergedBamFile $singlesBamFile $pairedBamFile");
        
        # Mark duplicates and use mpileup
        my $cleanedBam = $mappingDir . "/" . $readGroup . ".merged.cleaned.bam";
        my $sortedBam = $mappingDir . "/" . $readGroup . ".merged.cleaned.sorted.bam";
        my $markDupsBam = $mappingDir . "/" . $readGroup . ".merged.cleaned.sorted.markDups.bam";
        my $markDupsMetrics = $mappingDir . "/" . $readGroup . ".merged.cleaned.sorted.markDups.metrics";
        #my $pileupFile = $mappingDir . "/" . $readGroup . ".mpileup";
        system("java -Xmx16g -jar ~/bin/picard-tools-1.119/CleanSam.jar I=$mergedBamFile O=$cleanedBam");    
        system("java -Xmx16g -jar ~/bin/picard-tools-1.119/AddOrReplaceReadGroups.jar I=$cleanedBam O=$sortedBam SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=Lib1 RGID=$readGroup RGSM=$readGroup VALIDATION_STRINGENCY=LENIENT");
	system("java -Xmx16g -jar ~/bin/picard-tools-1.119/MarkDuplicates.jar I=$sortedBam O=$markDupsBam METRICS_FILE=$markDupsMetrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 ASSUME_SORTED=true REMOVE_DUPLICATES=false");
        # system("samtools mpileup $markDupsBam > $pileupFile");
    
        print "Stats for $readGroup:\n";
        system("samtools flagstat $markDupsBam");
        print "\n\n\n";
        # All we need to keep is the $markDupsBam file, so get rid of the other sams and bams
        unlink ($singlesBamFile, $pairedBamFile, $mergedBamFile, $cleanedBam); # Hold on to $sortedBam for now in case we don't actually want to mark duplicates
        
    }
}












#Documentation
__END__

=head1 NAME

tortReadQC.pl

=head1 SYNOPSIS 

perl tortReadQC.pl --reads <file> --adapters <file> --out <outputDirectory> --log <logfile.txt> --trim --map --threads 6 --reference <path/to/bwa/reference/genome/index.fasta>

 Options:
   -reads=s           Directory with raw reads in gzipped fastq format
   -adapters=s        Directory with adapters fasta files
   -out=s             Name of output directory (it will be created)
   -log=s             Name of logfile to print output (you will probably also want
                      to capture STDERR manually)
   -trim              (flag) Perform read trimming using Trimmomatic and join with fastq-join
   -map               (flag) Perform mapping to reference genome
   -threads=i         Threads to use for multithreading (default=4)
   -reference=s       Reference to use for mapping (default to /mnt/Data4/genomes/Galapagos.fasta)
   -help|man        Prints out documentation


=head1 DESCRIPTION

This was written for the purposes of QCing HiSeq reads for a desert tortoise
project

=cut