#! /usr/bin/env perl
#####Description of this code#####
#This is code is used to compare log10(Experimental_Overlap) values to the upper 99.9% one-sided confidence intervals (UP999).
###########################################################################################################
#####how to run this code #####
##command##
#BarBIQ_final_overlap_selete_by_PV_seq.pl --overlap file1 --PV file3 --EDrop file4
##explaination##
#file1: output file from BarBIQ_final_repeat_overlap_seq.pl
#file3: simulation_overlap, an 0.999 confidential line of the distribution of the simulated Poission overlap, the simulation was done by BarBIQ_add_Simulation_overlap_AB_Poisson.pl and BarBIQ_add_Simulation_overlap_up999.pl 
#file4: see EDrop.txt, ourput file from BarBIQ_final_fitting_OD.r
########################################################################################################## 
#####Install#####
##None
############################################################################################################

#####code##### 

use strict;
use warnings;
# print "Now you are running program: $0\n";
# print "Your parameters are: @ARGV\n";

my ($i,$overlapfile,$PV_file,$Edrop_file,$log_file);
# my $threshold=0.1; ## a threshold for deleting the bad dataset if the standard error of fitting is > 0.1 for this analysis
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--overlap")  {$overlapfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--PV") {$PV_file = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--EDrop") {$Edrop_file = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input: --overlap XXXX --PV XXXXXX --EDrop XXXXXXX (--log XXXXXXX)\"\n $!";}
    }
my $outputfile = "$overlapfile"."_PV_list"; ##
unlink $outputfile;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are running program: $0\n";
print "Your parameters are: @ARGV\n";

my %EDrop;
my @PV;
open(PV, $PV_file) || die "canot open input file '$PV_file' $!";
while(my $gi = <PV>)
   {
     chomp $gi;
     my @info=split(/\s+/,$gi);
     if($info[0] eq "A" && $info[1] eq "B") {
     print ("Your PV file $PV_file is OK!\n");
     }
    else{
     push @PV, [@info];
     }
   }
close PV;
open(ED, $Edrop_file) || die "canot open input file '$Edrop_file' $!";
my $SE;
my $EDrop;
my $sample_ID;
my $gi = <ED>;
chomp $gi;
my @info=split(/\s+/,$gi);
for (my $i=0; $i<=$#info; $i++)
   {
     if($info[$i] eq "ID") {$sample_ID = $i;}
     if($info[$i] eq "EDrop") {$EDrop = $i;}
     if($info[$i] eq "SE") {$SE = $i;}
   }
if(!($EDrop && $SE)) { die "your inputfile $Edrop_file is wrong!!!!0007\n"; }
while($gi = <ED>)
   {
     chomp $gi;
     my @info=split(/\s+/,$gi);
     if (exists $EDrop{$info[0]}) { die "your inputfile $Edrop_file is wrong!!!!\n"}
    else{
        # if($info[2] > $threshold) {print "$info[0] is a bad data for overlap analysis! The std err: $info[2]\n";}
        $EDrop{$info[0]} = 10**$info[1]; print "$info[0]\t$EDrop{$info[0]}\n";
        }
   }
close ED;

open(IN, $overlapfile) || die "canot open input file '$overlapfile' $!";
open(OUTF, '>>', $outputfile) or die "Could not open file '$outputfile' $!";

$gi = <IN>;
chomp $gi;
@info=split(/\s+/,$gi);
my $SEQ1;
for (my $i=0; $i<=$#info; $i++)
   {
     if($info[$i] eq "SEQ1") {$SEQ1 = $i;}
   }
my $repeats;
if($SEQ1)
   {
    $repeats=$SEQ1/5;
    print "You have $repeats repeats\n";
   }
else{die "your inputfile $overlapfile is wrong!!!!\n";}
my @sample_consider;
for (my $i=1; $i<=$repeats; $i++)
   {
     if($info[$i*5-5] =~ /S([0-9]+)_A/) 
        {  
          my $sample="S"."$1"; 
          my $id=$i*5-5;
          if(exists $EDrop{$sample}) { 
                push @sample_consider,[($id,$EDrop{$sample})];
                print OUTF ("$sample\t");
              }
         else{}
        }
     else{ die "your inputfile $overlapfile is wrong!!!!\n";}
   }
print OUTF ("SEQ1\tSEQ2\n");
while($gi = <IN>)
   {
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $keep=0;
     for (my $i=0; $i<=$#sample_consider; $i++)
         {
          my $overlap=$info[$sample_consider[$i][0]+4];
          my $x=$info[$sample_consider[$i][0]+1]*$info[$sample_consider[$i][0]+3]/$sample_consider[$i][1]*5000;
          my $threshold;
          for(my $j=0; $j<=$#PV; $j++)
             {
               if($x<=($PV[$j][0]*$PV[$j][1])) {$threshold = $PV[$j][3]; last;}
             }
          if($overlap>$threshold) {$keep=1;}
          if ($x == 0) {$keep="-";}
          if ($overlap == 1) {$keep=2;}
          print OUTF ("$keep\t");
          $keep=0;
         }
     print OUTF ("$info[$SEQ1]\t$info[$SEQ1+1]\n");
   }
   
close OUTF;
close IN;

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.19
