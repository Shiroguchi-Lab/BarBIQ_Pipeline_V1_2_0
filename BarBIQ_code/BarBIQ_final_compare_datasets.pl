#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used for compare the cell counts of cOTUs of all samples after removing the contaminated cOTUs.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_compare_datasets.pl sample_1 sample_2 ..... outputfile_name
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

##check the output name##
if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!! $!";}

##Get all different sequences for all samples##
my $head;
my %IDS;
for (my $i = 0; $i<=$#inputfile; $i++)
    {
     my $sample=$i+1;
     print "Compare sample $sample:\t$inputfile[$i]\n";  ##print the sample's input name
     open(FILE, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
     my $gi=<FILE>;
     chomp $gi;
     my @info=split(/\s+/,$gi);
     my $id;
     for(my $j=0; $j<= $#info; $j++)
         {
          if($info[$j] eq "ID") {$id=$j;}
          else{$head=$head."\t$info[$j]";}
         }
     if (! defined $id) {die "your inputfile $inputfile[$i] is wrong001!!!\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
               if(exists $IDS{$info[$id]}) {}
               else{
                    $IDS{$info[$id]}=$info[$id];
                   }  
         }
      close FILE; 
     }

## Merge datasets
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
my %out;
my %No;
foreach my $key (sort keys %IDS)
     {
      $out{$key}=$key;
     }
for (my $i = 0; $i<=$#inputfile; $i++)
     {
      open(FILE1, $inputfile[$i]) or die "cannot open input file '$inputfile[$i]' $!";
      undef %No;
      my $gi=<FILE1>;
      chomp $gi;
      my @info=split(/\s+/,$gi);
      my $samples=$#info;
      while($gi=<FILE1>)
           {
            chomp $gi;
            @info=split(/\s+/,$gi);
            my $IDID = shift @info;
            $No{$IDID} = join("\t", @info);
           }
      close FILE1;
      foreach my $key (sort keys %IDS)
           {
            if(exists $No{$key}) 
                {
                 $out{$key}="$out{$key}"."\t$No{$key}";
                }
            else{
                 for(my $s = 1; $s<=$samples; $s++)
                   {
                     $out{$key}="$out{$key}"."\t0";
                   }
                }
           }
      }
##print the head of the table
print OUTF ("ID$head\n");
##print the data and delete all zero data
foreach my $key (sort keys %out)
    {
      my @info=split(/\s+/,$out{$key});
      my $sum=0;
      for(my $s = 1; $s<=$#info; $s++)
          {
            $sum=$sum+$info[$s];
          }
      if($sum>0)
        {
         print OUTF ("$out{$key}\n");
        }
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
##2019.1.5

