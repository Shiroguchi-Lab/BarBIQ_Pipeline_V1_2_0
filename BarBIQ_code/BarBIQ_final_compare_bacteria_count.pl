#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used for compare the cell counts of cOTUs from different samples or replicates.
#Note that, for the next step, only the samples or replicates used the same blank control should be compared.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_compare_bacteria_count.pl sample_1 sample_2 ..... outputfile_name
##interpretation##
#sample_1: the input file name of sample 1
#sample_2: the input file name of sample 2 
#........: any other samples
#outputfile_name: a output file's name for saving the results
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

# my $repeats=$#inputfile+1;
# print "You have $repeats repeats in this experiment, so the suggested threshold to delete the possible junk sequencese are:\n";
 
##check the output name##
unlink $outputfile;
# if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!! $!";}
# my $outputfile_IF=$outputfile."_IF";
# unlink $outputfile_IF;

##Get all different sequences for all samples##
# open(OUTFIF, '>>', $outputfile_IF) or die "cannot open input file '$outputfile_IF' $!";
my ($gi, @info, %No, %ID, %seq, $i, @samples);
# print OUTFIF ("Index\tSample\tType\n");
my $samplename;
for ($i = 0; $i<=$#inputfile; $i++)
    {
     if ($inputfile[$i] =~ /_S([0-9]*)_/) {push @samples, "S$1"; $samplename=$`;}
     elsif ($inputfile[$i] =~ /\ABlank([0-9]*)/) {push @samples, "Blank$1"; $samplename="Blank$1";}
     else {die "Your input is wrong!!!\n";}
     print "Sammple_$samples[$i]:\t$inputfile[$i]\n";  ##print the sample's input name
     # print OUTFIF ("$samples[$i]\t$samplename\n");
     open(FILE, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
               if(exists $seq{$info[0]}) {}
               else{
                    $seq{$info[0]}=$info[0];
                   }  
         }
      close FILE; 
     }
# close OUTFIF;

##Renumber the sequences and print the number of eaah sequence for all samples and the mapping results to the outputfile 
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
my %out;
foreach my $key (sort keys %seq)
     {
      $out{$key}=$key;
     }
my $head;
for ($i = 0; $i<=$#inputfile; $i++)
     {
      # my $x=$i+1;
      # $head = $head."\tSample_$x"; 
      $head = $head."\t$samples[$i]";
     }
for ($i = 0; $i<=$#inputfile; $i++)
     {
      open(FILE1, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
      undef %No;
      undef %ID;
      while($gi=<FILE1>)
           {
            chomp $gi;
            @info=split(/\s+/,$gi);
            $No{$info[0]} = $info[1];
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
print OUTF ("ID$head\n");
foreach my $key (sort keys %out)
    {
      print OUTF ("$out{$key}\n");
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
##2018.8.8

