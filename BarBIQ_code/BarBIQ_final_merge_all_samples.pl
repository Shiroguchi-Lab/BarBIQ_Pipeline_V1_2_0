#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to merge all RepSeqs from all samples and to name the resulting unique RepSeq types BarBIQ-identified sequences (Bar sequences), and each was labelled by an ID number.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_merge_all_repeats_files.pl sample_1 sample_2 ..... outputfile_name
##explaination##
#sample_1: the input file name of sample 1
#sample_2: the input file name of sample 2
#........: other samples
#outputfile_name: output file's name for saving the results
####################################################################################################

#####code#######
use strict;
use warnings;
my @inputfile = @ARGV;
my $log_file = pop @inputfile;
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

##save the input file names and output file name##
# my @inputfile = @ARGV;
my $outputfile = pop @inputfile;

my $repeats=$#inputfile+1;
print "You totally have $repeats samples\n";
 
##check the output name##
if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!! $!";}

##Get all different sequences for all samples##
my ($gi, @info, %No, %ID, %seq, %link, $i);
for ($i = 0; $i<=$#inputfile; $i++)
    {
     my $x=$i+1;
     print "Sammple_$x:\t$inputfile[$i]\n";  ##print the sample's input name
     open(FILE, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
     $gi=<FILE>; chomp $gi; @info=split(/\s+/,$gi);
     my $s;
     my $ave;
     my $LINK;
     my $I1;
     my $R2;
     for (my $j = 0; $j<=$#info; $j++)
        {
           if($info[$j] eq "Sequence") { $s = $j; }
           if($info[$j] eq "AverageCount") { $ave = $j; }
           if($info[$j] eq "LINK") { $LINK = $j; }
           if($info[$j] eq "I1") { $I1 = $j; }
           if($info[$j] eq "R2") { $R2 = $j; }
        }
    if($s && $ave && $LINK && $I1 && $R2) {} else{die "error001\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
               if(exists $seq{$info[$s]}) {}
               else{
                    $seq{$info[$s]}=$info[$ave];
                    $link{$info[$s]}="$info[$LINK]\t$info[$I1]\t$info[$R2]";
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
      my $s;
      my $ave;
      $gi=<FILE1>; chomp $gi; @info=split(/\s+/,$gi);
      for (my $j = 0; $j<=$#info; $j++)
        {
           if($info[$j] eq "Sequence") { $s = $j; }
           if($info[$j] eq "AverageCount") { $ave = $j; }
        }
      while($gi=<FILE1>)
           {
            chomp $gi;
            @info=split(/\s+/,$gi);
            $No{$info[$s]} = $info[$ave];
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
print OUTF ("ID$head\tSequence\tLINK\tI1\tR2\n");
foreach my $key (sort keys %out)
    {
     print OUTF ("$out{$key}\t$key\t$link{$key}\n");
    }
close OUTF; 
print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.001
##2019.1.1

