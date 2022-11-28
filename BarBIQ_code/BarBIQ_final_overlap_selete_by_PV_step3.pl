#! /usr/bin/env perl
#####Description of this code#####
#This is code is used to selet the Bar-sequence pairs from the same bacterium.
################################################################################################################how to run this code #####
##command###
#BarBIQ_final_overlap_selete_by_PV_step3.pl --overlap file1 --select file2
##explaination##
#file1: output file from BarBIQ_final_repeat_overlap_seq.pl;
#file2: output file from BarBIQ_final_overlap_selete_by_PV_seq_threshold.pl
########################################################################################################## 
#####Install#####
##None
############################################################################################################
#####code##### 

use strict;
use warnings;
# print "now you are running program: $0\n";
# print "Your parameters are: @ARGV\n";

my ($i,$overlapfile,$sel_file,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--overlap")  {$overlapfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--select") {$sel_file = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input: --overlap XXXX --select XXXXXX (--log log_file)\"\n $!";}
    }


my $outputfile = "$overlapfile"."_select"; ##
unlink $outputfile;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "now you are running program: $0\n";
print "Your parameters are: @ARGV\n";

my %select;
open(SELECT, $sel_file) || die "canot open input file '$sel_file' $!";
my $gi = <SELECT>;
chomp $gi;
my $SEQ1;
my $Y;
my $N;
my @info=split(/\s+/,$gi);
for (my $i=0; $i<=$#info; $i++)
   {
     if($info[$i] eq "SEQ1") {$SEQ1 = $i;}
     if($info[$i] eq "Yes") {$Y = $i;}
     if($info[$i] eq "No") {$N = $i;}
   }
if($SEQ1 && $Y && $N) {}
else{die "your inputfile $sel_file is wrong!!!!\n";}
while($gi = <SELECT>)
   {
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $ratio=$info[$Y]/($info[$Y]+$info[$N]);
     my $pair = "$info[$SEQ1]_$info[$SEQ1+1]";
     $select{$pair} = "$info[$Y]\t$info[$N]\t$ratio";
     $pair = "$info[$SEQ1+1]_$info[$SEQ1]";
     $select{$pair} = "$info[$Y]\t$info[$N]\t$ratio";
   }
close SELECT;

open(IN, $overlapfile) || die "canot open input file '$overlapfile' $!";
open(OUTF, '>>', $outputfile) or die "Could not open file '$outputfile' $!";
$gi = <IN>;
chomp $gi;
@info=split(/\s+/,$gi);
undef $SEQ1;
for (my $i=0; $i<=$#info; $i++)
   {
     if($info[$i] eq "SEQ1") {$SEQ1 = $i;}
   }
if($SEQ1)
   {
   }
else{die "your inputfile $overlapfile is wrong!!!!\n";}
print OUTF ("$gi\tYes\tNo\tRatio\n");
while($gi = <IN>)
   {
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $pair = "$info[$SEQ1]_$info[$SEQ1+1]";
          if (exists $select{$pair}) {
          print OUTF ("$gi\t$select{$pair}\n");
         }
   }
   
close OUTF;
close IN;

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.20
