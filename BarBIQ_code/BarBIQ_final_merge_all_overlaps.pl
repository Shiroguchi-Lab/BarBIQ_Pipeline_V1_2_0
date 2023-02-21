#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to merge all overlaping data for identifying multiple Bar sequences from the same bacterium from different samples or replicates.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_repeat_overlap_seq.pl sample_1 sample_2 ..... outputfile_name
##explaination##
#sample_1: the input file name of sample 1
#sample_2: the input file name of sample 2
#........: other samples
#outputfile_name: output file's name for saving the results
####################################################################################################

#####code#######
use strict;
use warnings;


##save the input file names and output file name##
my @inputfile = @ARGV;

my $log_file = pop @inputfile;
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

my $outputfile = pop @inputfile;
# my $normalization = $outputfile."_normalization";
my $repeats =$#inputfile+1;
print "You have $repeats samples including replicates in this experiment\n";
##check the output name##
if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!! $!";}

my @samples;
##Get all different sequences for all samples##
my ($gi, @info, %No,%pair,%seq);
if ($inputfile[0] =~ /S([0-9]*)_/) {push @samples, "S$1"; } else {die "Your input is wrong!!!001\n";}
print "Sammple_$1:\t$inputfile[0]\n";
open(FILE, $inputfile[0]) or die "cannot open input file '$inputfile[0]' $!";

while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          my $p="$info[0]"."\t$info[2]";
               if(exists $pair{$p}) { die "Your input is wrong!!!005 $p\n";}
               else{
                    $pair{$p}="$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]";
                    $seq{$p} = "$info[5]\t$info[6]";
                   }
         }
for (my $i = 1; $i<=$#inputfile; $i++)
    {
     if ($inputfile[$i] =~ /S([0-9]*)_/) {push @samples, "S$1"; } else {die "Your input is wrong002!!!\n";}
     print "Sammple_$1:\t$inputfile[$i]\n";
     open(FILE, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          my $p="$info[0]"."\t$info[2]";
               if(exists $pair{$p}) { $pair{$p} = "$pair{$p}\t$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]";}
               else{
                    die "$p\tYour input is wrong!!!007\n";
                   }  
         }
      close FILE; 
     }

##Renumber the sequences and print the number of eaah sequence for all samples and the mapping results to the outputfile 
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
print OUTF ("$samples[0]_A_ID\t$samples[0]_A\t$samples[0]_B_ID\t$samples[0]_B\t$samples[0]_O");
for(my $i=1; $i<=$#samples; $i++)
     {
       print OUTF ("\t$samples[$i]_A_ID\t$samples[$i]_A\t$samples[$i]_B_ID\t$samples[$i]_B\t$samples[$i]_O");
     }
print OUTF ("\tSEQ1\tSEQ2\n");
foreach my $key (sort keys %pair)
     {
       print OUTF ("$pair{$key}\t$seq{$key}\n");
     }
close OUTF;
print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.19

