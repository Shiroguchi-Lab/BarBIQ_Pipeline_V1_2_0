#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
#This code is used to identify multiple Bar sequences (i.e., 16S rRNA sequences) for the same bacterium using the cellular barcodes by distinguishing a natural co-occurrence of Bar sequences of the same bacterium from an accidental co-occurrence of Bar sequences from different bacteria that existed in the same droplet.
#For this process, please use the processed data from BarBIQ_sub_clustering_step_two.pl using threshold 0.1.
##########################################################################################################################################
######how to run this code #####
###command##
##BarBIQ_final_droplets_overlap.pl --in inputfile --in inputfile
##explaination##
#in: the input file which is the output file from BarBIQ_final_review_seq_by_repeats.pl.
#repeat-file: the output file from BarBIQ_final_add_bacteria_name_to_repeats_LIB.pl. 
##########################################################################################################################################
#####Install#####
##None 
##########################################################################################################################################

#####code#####
use strict;
use warnings;
#use IPC::System::Simple qw(system);

#print "Now you are runing $0\n";
#print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$outputfile, $outputfile_cluster, $outputfile_IF, $repeats, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
#     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--bar") {$repeats = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --bar bar_sequence_file --log log_file\"\n $!";}
    }
if(!$log_file) {
                 $log_file = $inputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in inputfile\"\n $!";}
if(!$repeats)   {die "Your input is wrong!!!\n Please input \"--bar bar_sequence_file\"\n $!";}
# if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out: outputfile\"\n $!";}
$outputfile="$inputfile"."_DropOver";
# $outputfile_cluster="$inputfile"."_cluster";
# $outputfile_IF="$inputfile"."_IF";
unlink $outputfile;
# unlink $outputfile_cluster;
# unlink $outputfile_IF;
#if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfiles is:\n$inputfile\n";
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
if(!(($#info == 5) && ($info[0] =~ /\Acluster_/)))
   {
    die "Your input file $inputfile is wrong!!! please check!!!\n $!";
   }
close FILE;

if(!(-e $repeats)) {die "Your input file $repeats is not existed!!! please check!!!\n $!";}
open (FILE,$repeats) or die "Could not open file '$repeats' $!"; #
print "Your bar-sequence file is:\n$repeats\n";
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
my $LIBID;
my $seq;
for($i=0; $i<=$#info; $i++)
   {
    if($info[$i] eq "Bar-sequence-ID") { $LIBID = $i;}
    if($info[$i] eq "Sequence") { $seq = $i;}
   }
if (!(defined $LIBID && defined $seq)) {die "Your input file $repeats is wrong 002\n";}
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

###Main code###
my %LIBIDs;
while($gi=<FILE>)
    {
       chomp $gi;
       @info=split(/\s+/,$gi);
       $LIBIDs{$info[$seq]} = $info[$LIBID];
    }
close FILE;


my @sequence_cluster; ## statistic for each cluster
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
my ($cluster_name, $cluster_name_next);
if($info[0] =~ m{\Acluster_([0-9]*)})
        {
          $cluster_name = "cluster_"."$1";
        }
push @sequence_cluster, [@info];
my %cluster_seq; ## key is the cluster_name, valuse the lib IDs
my $IDIDID=":";
while($gi=<FILE>)
     {
       chomp $gi;
       @info=split(/\s+/,$gi);
       if($info[0] =~ m{\Acluster_([0-9]*)})
        {
          $cluster_name_next = "cluster_"."$1";
        }
       if($cluster_name_next eq $cluster_name)
           {
             push @sequence_cluster, [@info];
           }
       else{
            for(my $x=0; $x<=$#sequence_cluster; $x++)
                {
                 $IDIDID = "$IDIDID"."$LIBIDs{$sequence_cluster[$x][5]}".":";
                }
            $cluster_seq{$cluster_name}=$IDIDID;
            $cluster_name=$cluster_name_next;
            undef @sequence_cluster;
            $IDIDID=":";
            push @sequence_cluster, [@info];
           }
      }
for(my $x=0; $x<=$#sequence_cluster; $x++)
      {
       $IDIDID = "$IDIDID"."$LIBIDs{$sequence_cluster[$x][5]}".":";
      }
$cluster_seq{$cluster_name}=$IDIDID;
close FILE;

open (OUTF,'>>',$outputfile) or die "Could not open file '$outputfile' $!"; ## write the overlap results file
my @key = sort keys %LIBIDs;
my @overlapping;
for(my $i=0; $i<=$#key; $i++)
    {
      for(my $j=0; $j<=$#key; $j++)
          {
            my $seqi=":"."$LIBIDs{$key[$i]}".":";
            my $seqj=":"."$LIBIDs{$key[$j]}".":";
            my $No_two=0;
            foreach my $key (sort keys %cluster_seq)
               { 
                 if(($cluster_seq{$key} =~ /$seqi/) && ($cluster_seq{$key} =~ /$seqj/)) {$No_two++;}
               }
            $overlapping[$i][$j]=$No_two;
         }
     }

for (my $i=0; $i<$#overlapping; $i++)
   {
    for (my $j=$i+1; $j<=$#overlapping; $j++)
      {
        my $seqi=$LIBIDs{$key[$i]};
        my $seqj=$LIBIDs{$key[$j]};
            if ($seqi lt $seqj)
                {
                 print OUTF ("$seqi\t$overlapping[$i][$i]\t$seqj\t$overlapping[$j][$j]\t$overlapping[$i][$j]\t$key[$i]\t$key[$j]\n");
                }
          elsif ($seqj lt $seqi) {
                print OUTF ("$seqj\t$overlapping[$j][$j]\t$seqi\t$overlapping[$i][$i]\t$overlapping[$i][$j]\t$key[$j]\t$key[$i]\n");
               }
          else { die "Somthing is wrong008\n";}
      }
   }

close OUTF;
##statistic the information of clusters and mapping##

print "Done!!!\n";
###Main code###
####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.19
