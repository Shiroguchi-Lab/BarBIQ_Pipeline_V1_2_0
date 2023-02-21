#! /usr/bin/env perl
###############################################################
#######Description of this code#####
#This code is used to count the number of cells for each cOTU
##########################################################################
######how to run this code #####
###command###
#BarBIQ_final_count_bacteria.pl --in file1
###explaination###
#file1: output file from BarBIQ_final_add_bacteria_name_to_final_rep_seq_file.pl
############################################
#####Install#####
#None
###################################################

#####code#####

use strict;
use warnings;

# print "Now you are runing \n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$outputfile,$outputfile_cluster,$outputfile_count, $outputfile_IF,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
      elsif($ARGV[$i] eq "--log")  {$log_file = $ARGV[$i+1];}
    # elsif ($ARGV[$i] eq "--group") {$bacteria_name = $ARGV[$i+1];}
#     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];} 
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --log log_file\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
# if(!$bacteria_name)   {die "Your input is wrong!!!\n Please input \"--name bacteria_name_file\"\n $!";}
# if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out: outputfile\"\n $!";}
$outputfile_count="$inputfile"."_bac_count";
unlink $outputfile_count;

if(!$log_file) {
                 $log_file = $outputfile_count."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

# $outputfile_IF="$inputfile"."_bac_IF";
$outputfile_cluster="$inputfile"."_bac_cell";
unlink $outputfile_cluster;
if(-e $outputfile_count){die "Your output file $outputfile_count is already existed!!! please check!!!\n $!";}
# if(-e $outputfile_IF){die "Your output file $outputfile_IF is already existed!!! please check!!!\n $!";}
if(-e $outputfile_cluster){die "Your output file $outputfile_cluster is already existed!!! please check!!!\n $!";}
##read command##
#
##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfiles is:\n$inputfile\n";
my $gi=<FILE>; chomp $gi;
my @info=split(/\s+/,$gi);
if(!(($info[6] =~ /\ADroplet_/) && ($info[0] =~ /\Acluster_/)))
   {
    die "Your input file $inputfile is wrong!!! please check!!!\n $!";
   }
close FILE;
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##
#
###Main code###

open(FILE, $inputfile) or die "cannot open input file '$inputfile' $!";
my %Bacteria;
my @clusters;
while($gi=<FILE>)
{
     chomp $gi;
     @info=split(/\s+/,$gi);
     push @clusters, [($info[6], $info[8], $info[1])];
     $Bacteria{$info[8]}=0;
}
close FILE;

my @clusters_sort = sort { $a->[0] cmp $b->[0]
                       or $a->[1] cmp $b->[1]}
                       @clusters;
my %clusters;
my $reads;
my $cluster_ID=$clusters_sort[0][0];
my $bacteria_ID=$clusters_sort[0][1];
foreach (keys %Bacteria)
    {
     $Bacteria{$_}=0;
    }
$Bacteria{$bacteria_ID}++;
$clusters{$cluster_ID} = ":"."$bacteria_ID".":";
open(OUTF1, '>>', $outputfile_cluster)  or die "canot open input file '$outputfile_cluster' $!";
print OUTF1 ("Droplet-ID\tcOTU-ID\tNo-of-reads\n");
$reads=$clusters_sort[0][2];
for ($i=1; $i<=$#clusters_sort; $i++)
    {
     if($clusters_sort[$i][0] eq $cluster_ID)
         {
          if($clusters_sort[$i][1] eq $bacteria_ID)
              {
               $reads = $reads+$clusters_sort[$i][2];
              }
          else{
               print OUTF1 ("$cluster_ID\t$bacteria_ID\t$reads\n");
               $bacteria_ID=$clusters_sort[$i][1];
               $Bacteria{$bacteria_ID}++;
               $reads=$clusters_sort[$i][2];
               $clusters{$cluster_ID} = "$clusters{$cluster_ID}".":$bacteria_ID".":";
              }
         }
     else{
          print OUTF1 ("$cluster_ID\t$bacteria_ID\t$reads\n");
          $cluster_ID=$clusters_sort[$i][0];
          $bacteria_ID=$clusters_sort[$i][1];
          $reads=$clusters_sort[$i][2];
          $Bacteria{$bacteria_ID}++;
          $clusters{$cluster_ID} = ":"."$bacteria_ID".":";
         }
      }
print OUTF1 ("$cluster_ID\t$bacteria_ID\t$reads\n");
close OUTF1;

open(OUTF2, '>>', $outputfile_count)  or die "canot open input file '$outputfile_count' $!";
foreach (keys %Bacteria)
    {
     print OUTF2 ("$_\t$Bacteria{$_}\n");
    }
close OUTF2;

# open (IF,'>>',$outputfile_IF) or die "Could not open file '$outputfile_IF' $!";
# my @key = keys %Bacteria;
# my @overlapping;
#for(my $i=0; $i<=$#key; $i++)
#    {
#      for(my $j=0; $j<=$#key; $j++)
#          { my $seqi=":"."$key[$i]".":";
#            my $seqj=":"."$key[$j]".":";
#            my $No_two=0;            
#            foreach my $key (sort keys %clusters)
#               {
#                 if(($clusters{$key} =~ /$seqi/) && ($clusters{$key} =~ /$seqj/)) {$No_two++;}
#               }
#            $overlapping[$i][$j]=$No_two;
#          }
#     }
#my $out;
#my $string=join("\t",@key);
#print IF ("$string\n");
#for $out(@overlapping)
#  {
#   $string=join("\t",@$out);
#   print IF ("$string\n");
#  }
#close IF;

print "Done!!!\n";
###Main code###
#####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.20




