#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to remove low-count RepSeq types because unexpected errors might occur.
#Different thresholds are applied for different sampling replicates.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_repeat_delete_low_count.pl inputfile 
##interpretation##
#inputfile: the output file from BarBIQ_final_repeat_delete_low_count.pl
####################################################################################################

#####code#######
use strict;
use warnings;

my $log_file = $ARGV[1];
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

##save the input file names and output file name##
my $inputfile = $ARGV[0];
my $outputfile = $inputfile."_del_low";
unlink $outputfile;
open(FILE, $inputfile) or die "cannot open input file '$inputfile' $!";
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
##print the head of the table
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
my $ave;
my $seq;
for(my $i=1; $i<=$#info; $i++)
   {
    if($info[$i] eq "AverageCount") {$ave=$i;}
    if($info[$i] eq "Sequence") {$seq=$i;}
   }
if (!($ave && $seq)) {die "Your input file $inputfile is wrong!!! please check!!!\n $!";}
my $repeats=$seq-1;
print "You have $repeats repeats in this experiment, so the suggested threshold to delete the possible junk sequencese are:\n";
my $delete_threshold=6/$repeats; ## we suggest that the  threshold to delete the possible junk sequencese should be 6/No. of repeats. When your repeats are <=3.
# $delete_threshold=10; ##you can set your sepecific threshold here.  
print ">=$delete_threshold\n";
print OUTF ("$gi\n");
while($gi=<FILE>)
  {
     chomp $gi;
     @info=split(/\s+/,$gi);
     if($info[$ave] >= $delete_threshold)
      {
       print OUTF ("$gi\n");
      }
  }        

close OUTF; 
close FILE;


print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.18

