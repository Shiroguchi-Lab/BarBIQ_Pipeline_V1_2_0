#! /usr/bin/env perl
#####Description of this code#####
#This is code is used to identify Bar-sequence pairs from the same bacterium, by a threshold of Ratio_Positive > 0.5.
################################################################################################################how to run this code #####
##command###
#BarBIQ_final_overlap_selete_by_PV_seq_threshold.pl inputfile
##explaination##
#inputfile: output file from BarBIQ_final_overlap_selete_by_PV_seq.pl
########################################################################################################## 
#####Install#####
##None
############################################################################################################
#####code##### 

use strict;
use warnings;

my $log_file = $ARGV[2];
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are running program: $0\n";
print "Your parameters are: @ARGV\n";

my $inputfile="$ARGV[0]";
my $outputfile = "$inputfile"."_$ARGV[1]"; ##
unlink $outputfile;
my $threshold=$ARGV[1];
open(FILE, $inputfile) || die "canot open input file '$inputfile' $!";
open(OUTF, '>>', $outputfile) or die "Could not open file '$outputfile' $!";
my $gi = <FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
my $SEQ1;
for (my $i=0; $i<=$#info; $i++)
   {
     if($info[$i] eq "SEQ1") {$SEQ1 = $i;}
   }
if($SEQ1)
   {
    print "You have $SEQ1 repeats\n";
   }
else{die "your inputfile $inputfile is wrong!!!!\n";}
print OUTF ("$gi\tYes\tNo\n");
while($gi = <FILE>)
   {
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $no=0;
     my $yes=0;
     for (my $i=0; $i<$SEQ1; $i++)
         {
          if($info[$i] eq "0") {$no++;}
          if($info[$i] eq "1") {$yes++;}
         }
   if(($no+$yes)>1)
     {
     if($yes/($no+$yes) > $threshold) {
             print OUTF ("$gi\t$yes\t$no\n");
            }
     }
   }
   
close OUTF;
close FILE;
##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.19
