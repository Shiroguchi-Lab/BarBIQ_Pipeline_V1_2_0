#! /usr/bin/env perl
###############################################################
#######Description of this code#####
#This code is used to add the cOTU IDs to the BCluser and RepSeq file from BarBIQ_final_review_seq_by_repeats.pl.
##########################################################################
######how to run this code #####
###command###
#BarBIQ_final_add_bacteria_name_to_final_rep_seq_file.pl --in file1 --COTU file2
###explaination###
#file1: output file from BarBIQ_final_review_seq_by_repeats.pl
#file2: output file from BarBIQ_final_lib_COTU_ID.pl
############################################
#####Install#####
#None
###################################################

#####code#####

use strict;
use warnings;

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$lib_COUT,$inputfile, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--COTU") {$lib_COUT = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --COTU lib_COTU (--log log_file)\"\n $!";}
    }
if(!$lib_COUT)   {die "Your input is wrong!!!\n Please input \"--repeat repeatfile\"\n $!";}
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
my $outputfile = $inputfile."_ID";
my $outputfile_stat=$outputfile."_stat";
unlink $outputfile;
unlink $outputfile_stat;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

##read command##
#
##check the inputfile##
if(!(-e $lib_COUT)) {die "Your input file $lib_COUT is not existed!!! please check!!!\n $!";}
open (FILE,$lib_COUT) or die "Could not open file '$lib_COUT' $!"; # open inputfile
print "Your lib_COUT file is:\n$lib_COUT\n";
my $gi=<FILE>; chomp $gi;
my @info=split(/\s+/,$gi);
if(!($info[1] eq "cOTU-ID"))
   {
    die "Your input file $lib_COUT is wrong!!! please check!!!001\n $!";
   }
$gi=<FILE>; chomp $gi;
@info=split(/\s+/,$gi);
if(!($info[1] =~ /\ACOTU/))
   {
    print $info[1];
    die "Your input file $lib_COUT is wrong!!! please check005!!!\n $!";
   }
close FILE;
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##
#
###Main code###
open(FILE, $lib_COUT) or die "cannot open input file '$lib_COUT' $!";
$gi=<FILE>;chomp $gi;
my $COUT_ID;
my $seq;
my %ID_bac;
my %ID_seq;
my $seq_ID;
my %stat;
@info=split(/\s+/,$gi);
for($i=0; $i<=$#info; $i++)
    {
      if ($info[$i] eq "Bar-sequence-ID") { $seq_ID=$i;}
      if ($info[$i] eq "cOTU-ID") { $COUT_ID=$i;}
      if ($info[$i] eq "Sequence") { $seq=$i;}
    }
if (!(defined $seq_ID && defined $seq && defined $COUT_ID)) {die "Your input file $lib_COUT is wrong!!!\n";}
while($gi=<FILE>)
   {
     chomp $gi;
     @info=split(/\s+/,$gi);
     $ID_bac{$info[$seq]}=$info[$COUT_ID];
     $ID_seq{$info[$seq]}=$info[$seq_ID];
     $stat{$info[$seq_ID]} = 0;
 #     print "$info[$seq]\n$info[$name]\n";
   }
close FILE;


open(FILE, $inputfile) or die "cannot open input file '$inputfile' $!";
open(OUTF, '>>', $outputfile)  or die "canot open input file '$outputfile' $!";
while($gi=<FILE>)
{
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $cluster_name;
     if($info[0] =~ m{\Acluster_([0-9]*)})
        {
          $cluster_name = "Droplet_"."$1";
        }
     if (exists $ID_bac{$info[5]})
        {
#         print "$info[5]\n";
         print OUTF ("$gi\t$cluster_name\t$ID_seq{$info[5]}\t$ID_bac{$info[5]}\n");
         $stat{$ID_seq{$info[5]}}++;
        }
    else{
         print "$info[5]\n";
         die "Your input file $inputfile is wrong!!!\n";
        }
}
close FILE;
close OUTF;

open(OUTF, '>>', $outputfile_stat)  or die "canot open input file '$outputfile_stat' $!";
print OUTF ("Bar-sequence-ID\tNo-of-reads\n");
foreach my $key (sort keys %stat)
  {
    print OUTF ("$key\t$stat{$key}\n");
  }

close OUTF;

print "Done\n";
###Main code###
#####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.20




