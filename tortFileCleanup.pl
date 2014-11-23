#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my $help = 0;
my $threads = 32;


GetOptions  ("help|man"     => \$help,
             "threads=i"    => \$threads) || pod2usage(2);

if ($help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

my @samples = ("etort-1","etort-2","etort-3","etort-4","etort-5","etort-6","etort-7","etort-8","etort-9","etort-10","etort-11","etort-12","etort-13","etort-14","etort-15","etort-16","etort-17","etort-18","etort-19","etort-20","etort-21","etort-22","etort-23","etort-24","etort-25","etort-26","etort-27","etort-28","etort-30","etort-31","etort-32","etort-33","etort-34","etort-35","etort-36","etort-37","etort-38","etort-39","etort-40","etort-41","etort-42","etort-43","etort-44","etort-45","etort-46","etort-47","etort-48","etort-49","etort-50","etort-51","etort-52","etort-53","etort-54","etort-55","etort-56","etort-57","etort-58","etort-59","etort-60","etort-61","etort-62","etort-63","etort-64","etort-65","etort-66","etort-67","etort-68","etort-69","etort-70","etort-71","etort-72","etort-73","etort-74","etort-75","etort-76","etort-77","etort-78","etort-79","etort-80","etort-81","etort-82","etort-83","etort-84","etort-85","etort-86","etort-87","etort-88","etort-89","etort-90","etort-91","etort-92","etort-93","etort-94","etort-95","etort-96","etort-97","etort-98","etort-99","etort-100","etort-102","etort-103","etort-104","etort-105","etort-106","etort-107","etort-108","etort-109","etort-110","etort-111","etort-112","etort-113","etort-115","etort-116","etort-117","etort-119","etort-120","etort-121","etort-122","etort-123","etort-124","etort-125","etort-126","etort-127","etort-128","etort-129","etort-130","etort-131","etort-132","etort-133","etort-134","etort-135","etort-136","etort-137","etort-138","etort-139","etort-140","etort-141","etort-142","etort-143","etort-144","etort-145","etort-146","etort-147","etort-148","etort-149","etort-150","etort-151","etort-152","etort-153","etort-154","etort-155","etort-156","etort-157","etort-158","etort-159","etort-161","etort-162","etort-163","etort-164","etort-166","etort-167","etort-168","etort-169","etort-170","etort-172","etort-173","etort-174","etort-175","etort-176","etort-177","etort-178","etort-179","etort-180","etort-181","etort-182","etort-183","etort-184","etort-185","etort-186","etort-187","etort-188","etort-189","etort-190","etort-191","etort-192","etort-193","etort-194","etort-195","etort-196","etort-197","etort-202","etort-203","etort-208","etort-209","etort-210","etort-211","etort-212","etort-213","etort-214","etort-215","etort-216","etort-217","etort-218","etort-219","etort-220","etort-221","etort-222","etort-223","etort-224","etort-225","etort-226","etort-227","etort-228","etort-229","etort-230","etort-231","etort-232","etort-233","etort-234","etort-235","etort-236","etort-237","etort-238","etort-239","etort-240","etort-241","etort-242","etort-243","etort-244","etort-245","etort-247","etort-248","etort-249","etort-250","etort-251","etort-252","etort-253","etort-254","etort-255","etort-256","etort-257","etort-258","etort-259","etort-260","etort-261","etort-262","etort-263","etort-264","etort-265","etort-266","etort-267","etort-268","etort-269","etort-270","etort-271","etort-272","etort-273","etort-274","etort-275","etort-276","etort-277","etort-278","etort-282","etort-283","etort-284","etort-285","etort-286","etort-287","etort-288","etort-289","etort-296","etort-297","etort-298","etort-299","etort-300","etort-301");

# Remove both singles files and the joined files (because we have the concatenated version as one single-end file)
my $unlinkForkManager = Parallel::ForkManager($threads);
foreach my $tort (@samples) {
    $unlinkForkManager->start and next;
    my $singles1File = "qc/trimmomatic/$tort" . "_R1s_Ns_trim.fastq.gz";
    my $singles2File = "qc/trimmomatic/$tort" . "_R2s_Ns_trim.fastq.gz";
    my $joinedFile = "qc/fastq-join/$tort" . "_Ns_trim_fqj.join.fastq.gz";
    unlink($singles1File);
    unlink($singles2File);
    unlink($joinedFile);
    $unlinkForkManager->finish;
}
$unlinkForkManager->wait_all_children;






#Documentation
__END__

=head1 NAME

tortFileCleanup.pl

=head1 SYNOPSIS 

perl tortFileCleanup.pl

 Options:
   -threads=i       Number of threads to use
   -help|man        Prints out documentation


=head1 DESCRIPTION

This script removes redundant read files (singletons and joined files) for all samples

=cut