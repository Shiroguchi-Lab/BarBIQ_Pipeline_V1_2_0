#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to calculate averaged counts for each RepSeq from different measurements (sampling replicates).
#To calculate the average counts, the counts of RepSeq types for each replicate were normalized by the total count of all RepSeq types, and the total count after normalization was the same as the highest unnormalized total count among three replicates. 
#the inputfiles are the output file from BarBIQ_M_I1R2.pl or BarBIQ_sub_statistic.pl.
#This code can calculate from any No. of replicates.
#The output file includes renumbered RepSeq ID, count in each replicate, average counts from all replicates for RepSeqs.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_repeats_no_junk_deletion.pl repeat_1 repeat_2 ..... outputfile_name
##explaination##
#repeat_1: the input file name of replicate 1 of a sample(this code can only calculate one sample)
#repeat_2: the input file name of replicate 2 of a sample
#........: other replicates
#outputfile_name: output file's name for saving the results
####################################################################################################

#####code#######
use strict;
use warnings;

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##save the input file names and output file name##
my @inputfile = @ARGV;
my $log_file = pop @inputfile;
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

my $outputfile = pop @inputfile;
my $normalization = $outputfile."_normalization";
unlink $normalization;

my $repeats=$#inputfile+1;
my $delete_threshold=0; ## meaning do not delete any RepSeqs
 

##Get all different sequences for all samples##
my ($gi, @info, %No, %ID, %seq, $i);
for ($i = 0; $i<=$#inputfile; $i++)
    {
     my $x=$i+1;
     print "Sammple_$x:\t$inputfile[$i]\n";  ##print the sample's input name
     open(FILE, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          if($info[3] eq "LINK" || $info[3] eq "I1R2") {}
        else{die "Your inputfile is wrong!!!003\n";}
          if(exists $seq{$info[2]}) {}
        else{
             $seq{$info[2]}="$info[3]\t$info[4]\t$info[5]";
            }  
         }
      close FILE; 
     }

##Renumber the sequences and print the number of eaah sequence for all samples and the mapping results to the outputfile 
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
my $seq=100000;
my %out;
foreach my $key (sort keys %seq)
     {
      $seq++;
      $out{$key}="sequence_$seq";
     }
my $head;
for ($i = 0; $i<=$#inputfile; $i++)
     {
      # my $x=$i+1;
      # $head = $head."\tSample_$x"; 
      $head = $head."\t$inputfile[$i]";
     }
for ($i = 0; $i<=$#inputfile; $i++)
     {
      open(FILE1, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
      undef %No;
      while($gi=<FILE1>)
           {
            chomp $gi;
            @info=split(/\s+/,$gi);
            $No{$info[2]} = $info[1];
           }
      close FILE1;
      foreach my $key (sort keys %seq)
           {
            if(exists $No{$key}) 
                {
                 $out{$key}="$out{$key}"."\t$No{$key}";
                }
            else{
                 $out{$key}="$out{$key}"."\t0";
                }
           }
      }
##print the head of the table
print OUTF ("ID$head\tSequence\tLINK\tI1\tR2\tAverageCount\n");
my @sum;
for ($i = 0; $i<=$#inputfile; $i++)
   {
      $sum[$i+1]=0;
   }
foreach my $key (sort keys %out)
    {
     my @info=split(/\s+/,$out{$key});
     my $sum=0;
     my $zero=1;
     for(my $i=1; $i<=$#info; $i++)
        {
          $sum=$sum+$info[$i];
          $zero=$zero*$info[$i];
        }
     my $ave=$sum/$repeats;
     if ($zero != 0) 
         {
            for(my $i=1; $i<=$#info; $i++)
               {
                  $sum[$i]=$sum[$i]+$info[$i];
               }
         }
     if ($ave>=$delete_threshold)
            {
             print OUTF ("$out{$key}\t$key\t$seq{$key}\t$ave\n");
            }
    }
close OUTF; 
my $max_sum=0;
for ($i = 0; $i<=$#inputfile; $i++)
   {
      print "$sum[$i+1]\n";
      if($max_sum<$sum[$i+1]){$max_sum=$sum[$i+1];}
   }
print "The max:$max_sum\n";
open(NORM, '>>', $normalization) or die "cannot open input file '$normalization' $!";
print NORM ("ID$head\tSequence\tLINK\tI1\tR2\tAverageCount\n");

foreach my $key (sort keys %out)
    {
     my @info=split(/\s+/,$out{$key});
     my $sum=0;
     for(my $i=1; $i<=$#info; $i++)
        {
          $info[$i] = $info[$i]/$sum[$i]*$max_sum;
          $sum=$sum+$info[$i];
        }
     my $ave=$sum/$repeats;
     my $string=join("\t",@info);
     if ($ave>=$delete_threshold)
            {
             print NORM ("$string\t$key\t$seq{$key}\t$ave\n");
            }
    }
close NORM;
unlink $outputfile;

print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.18

