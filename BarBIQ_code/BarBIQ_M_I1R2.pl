#! /usr/bin/env perl
##########################################################################################################
#####Description of this code#####
# This code is used to cluster the 16S rRNA sequences depends on I1 and R2 reads within each BCluster (clustered by the barcode, R1 reads), and to get the representative sequences (RepSeqs) for each obtained sub-cluster(SCluster);
# and also to remove shifted RepSeqs;
# and also to link I1 and R2 RepSeqs;
# and also to remove one insertion and deletion (1-Indel) RepSeqs;
# and also to remove chimaeras;
# and also to remove other Indel-related errored RepSeqs;
# finally, to count BClusters for each RepSeq in each index.
#### As an option, you can only process I1 or R2 reads alone.
# It is a combined pipline of 
#    BarBIQ_sub_clean_end_primer_and_barcode.pl 
#    BarBIQ_sub_clustering_step_one.pl
#    BarBIQ_sub_clustering_step_two.pl 
#    BarBIQ_sub_shift.pl
#    BarBIQ_sub_link.pl
#    BarBIQ_sub_link_error.pl
#    BarBIQ_sub_indels.pl
#    BarBIQ_sub_chimeras.pl
#    BarBIQ_sub_similar.pl
#    BarBIQ_sub_statistic.pl
# So please make sure these codes are all under the same directory of this combined code.
# If you using this combined code instead of those separated codes above respectively, meaning you want to use all default parameters without any modification.
# The default parameters please check each code in detail. 
# The input data should be fastq or fastq.gz formated files from Miseq directly or unzipped. 
# Output file will be the sequences of the RepSeqs and the counts of these RepSeqs.
###########################################################################################################
#####how to run this code #####
##command##
# BarBIQ_M_I1R2.pl --I1 I1.fastq --R2 R2.fastq --out outputfile --bar R1file_bar --end-I1 No --end-R2 No. (--primer-I1 No. --primer-R2 No. --single Yes/No --middle Yes/No)
##explanation##
# --I1: I1 reads file, which is obtained directly from illumina Miseq machine, should be fastq format, end with .fastq or .fastq.gz
# --R2: R2 reads file, which is obtained directly from illumina Miseq machine, should be fastq format, end with .fastq or .fastq.gz * The order of all reads in I1 and R2 files should be the same
# --out: outputfile, please set a name for your outputfile, which will be used to generate the processing middle file names and the final output file name.
# --bar: Clustered R1 reads file, which is the outputfile from BarBIQ_index_and_leakage.pl or BarBIQ_M_R1.pl, please check this R1 file is from the same sample of the I1&R2 reads files
# --end-I1: a position number from the 5' end for trimming the 3'end of the I1 reads, you may decide it using code BarBIQ_sub_average_quality_of_each_position.pl or BarBIQ_M_R1.pl
# --end-R2: a position number from the 5' end for trimming the 3'end of the R2 reads, you may decide it using code BarBIQ_sub_average_quality_of_each_position.pl or BarBIQ_M_R1.pl
# --primer-I1: the length of the primer which used for amplifying the 16S rRNA sequences at the I1 read side, in our case we used 21 bases, so the default is 21, but you can change using --primer-I1 No.
# --primer-R2: the length of the primer which used for amplifying the 16S rRNA sequences at the R2 read side, in our case we used 17 bases, so the default is 17, but you can change using --primer-R2 No.
# --delMF: you may delete the middle files gerenated during processing using this parameter. defaut is "No" 
########################################################################################################## 
#####Install#####
## Please install the nucleotide-sequence-clusterizer and add it to Environment variable(let it possibale to be called directly)
## please install the perl Modules, IPC::System::Simple, and Text::Levenshtein::XS, before using this code
############################################################################################################

use strict;
use warnings;
use IPC::System::Simple qw(system);
# print "Now you are running program: $0\n";
# print "Your parameters are: @ARGV\n";
# print "Started at: ";
# print scalar localtime;
# print "\n";

##read command##
my $keep_middle="No";
my ($i,$inputfile_I1,$inputfile_R2,$outputfile,$cluster_R1,$end_I1,$end_R2,$pippath,$software,$log_file);
# my $single="No"; # default will process pair end data
my $primer_I1=21; # default primer length at the I1 side
my $primer_R2=17; # default primer length at the R2 side
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--I1") {$inputfile_I1 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--R2") {$inputfile_R2 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--bar") {$cluster_R1 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--end-I1") {$end_I1 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--end-R2") {$end_R2 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--primer-I1") {$primer_I1 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--primer-R2") {$primer_R2 = $ARGV[$i+1];}
  #   elsif ($ARGV[$i] eq "--single") {$single = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--middle") {$keep_middle = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--pip") {$pippath = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--soft") {$software = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else  {die "Your input is wrong!!!\n Please input: \"BarBIQ_M_I1R2.pl --I1 I1.fastq --R2 R2.fastq --out outputfile --bar R1file_bar --end-I1 NO --end-R2 NO --pip pipline-path --log log_file_name (option: --primer-I1 No --primer-R2 No --middle Yes/No --log log_file) --soft nucleotide-sequence-clusterizer\"\n $!";}
    }
if(!$log_file)  {die "Your input is wrong!!!\n Please input \"--log log_file_name\"\n $!";}
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are running program: $0\n";
print "Your parameters are: @ARGV\n";
print "Started at: ";
print scalar localtime;
print "\n";

if(!$inputfile_I1)   {die "Your input is wrong!!!\n Please input: \"--I1: I1.fastq\"\n $!";}
if ($inputfile_I1 =~ m{.fastq.gz\Z})
    {
        print "$inputfile_I1 is unzipping...\n";
        my $inputfile_I1_unzip=$`.".fastq";
        system "gunzip -c $inputfile_I1 > $inputfile_I1_unzip"; print "unzip finished!!!\n";
        $inputfile_I1=$`.".fastq";
    }
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out outputfile\"\n $!";}
#if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
if(!$cluster_R1)   {die "Your input is wrong!!!\n Please input \"--bar R1file_cluster\"\n $!";}
if(!$end_I1)   {die "Your input is wrong!!!\n Please input \"--end-I1 NO\"\n $!";}
      if(!$inputfile_R2)   {die "Your input is wrong!!!\n Please input \"--R2: R2.fastq\"\n $!";}
      if(!$end_R2)   {die "Your input is wrong!!!\n Please input \"--end-R2 NO\"\n $!";}
      if ($inputfile_R2 =~ m{.fastq.gz\Z})
         {
          print "$inputfile_R2 is unzipping...\n";
          my $inputfile_R2_unzip=$`.".fastq";
          system "gunzip -c $inputfile_R2 > $inputfile_R2_unzip"; print "unzip finished!!!\n";
          $inputfile_R2=$`.".fastq";
         }
if(!$pippath) {die "Your input is wrong!!!\n Please input \"--pip XXXX\"\n $!";}
##Finished read command##

##check the inputfile##
if(!($inputfile_I1 =~ m{.fastq\Z})){die "Your inputfile '$inputfile_I1' is not fastq!!! please check!!! $!";}
open (FILE,$inputfile_I1) or die "Could not open file '$inputfile_I1' $!"; # open inputfile
my $gi;
chomp ($gi=<FILE>);
   if($gi =~ m/\A@(.*)/s) # check the first line
       {
        $gi=<FILE>;$gi=<FILE>;$gi=<FILE>;chomp ($gi=<FILE>); # read to the 5th line
        if(!($gi =~ m/\A@(.*)/s)) # check the 5th line
            {die "Your inputfile '$inputfile_I1' is not fastq!!! please check!!! $!";}
       }
 else  {die "Your inputfile '$inputfile_I1' is not fastq format!!! please check!!! $!";}
close FILE; # close inputfile
     if(!($inputfile_R2 =~ m{.fastq\Z})){die "Your inputfile '$inputfile_R2' is not fastq!!! please check!!! $!";}
     open (FILE,$inputfile_R2) or die "Could not open file '$inputfile_R2' $!"; # open inputfile
     chomp ($gi=<FILE>);
     if($gi =~ m/\A@(.*)/s) # check the first line
         {
          $gi=<FILE>;$gi=<FILE>;$gi=<FILE>;chomp ($gi=<FILE>); # read to the 5th line
          if(!($gi =~ m/\A@(.*)/s)) # check the 5th line
              {die "Your inputfile '$inputfile_R2' is not fastq!!! please check!!! $!";}
         }
     else  {die "Your inputfile '$inputfile_R2' is not fastq format!!! please check!!! $!";}
     close FILE; # close inputfile
open (FILE,$cluster_R1) or die "Could not open file '$cluster_R1' $!"; # open inputfile
chomp ($gi=<FILE>);
my @info=split(/\s+/,$gi);
if(!($info[3] =~ m{\Acluster_})){die "Your inputfile '$cluster_R1' is wrong!!! please check!!! $!";}
print "Inputfile is OK!\nStart to calculating:\n";
##Finish check the inputfile##

##Run BarBIQ_sub_clean_clean_end_primer_and_barcode.pl##
my $outputfile_BarBIQ1="$outputfile"."_clean";
unlink $outputfile_BarBIQ1;
my $code = $pippath.'BarBIQ_sub_clean_end_primer_and_barcode.pl';
   system $code, '--I1',$inputfile_I1, '--R2', $inputfile_R2, '--out', $outputfile_BarBIQ1, '--bar', $cluster_R1, '--end-I1',  $end_I1, '--end-R2', $end_R2, '--primer-I1', $primer_I1, '--primer-R2', $primer_R2, "--log", $log_file;
##Finished BarBIQ_sub_clean_clean_end_primer_and_barcode.pl##

##Run BarBIQ_sub_clustering_step_one.pl##
my $outputfile_BarBIQ2="$outputfile_BarBIQ1"."_sub1";
unlink $outputfile_BarBIQ2;
$code = $pippath.'BarBIQ_sub_clustering_step_one.pl';
system $code,'--in',$outputfile_BarBIQ1,'--out',$outputfile_BarBIQ2, '--soft', $software, '--log', $log_file;
##Finished BarBIQ_sub_clustering_step_one.pl##

##Run BarBIQ_sub_clustering_step_two.pl##
my $outputfile_BarBIQ3="$outputfile_BarBIQ2"."_sub2";
unlink $outputfile_BarBIQ3;
$code = $pippath.'BarBIQ_sub_clustering_step_two.pl';
system $code,'--in',$outputfile_BarBIQ2,'--out', $outputfile_BarBIQ3, '--log', $log_file;
##Finished BarBIQ_sub_clustering_step_two.pl##

##Run BarBIQ_sub_shift.pl##
my $outputfile_BarBIQ4="$outputfile_BarBIQ3"."_shift";
unlink $outputfile_BarBIQ4;
$code = $pippath.'BarBIQ_sub_shift.pl';
   system $code,'--in', $outputfile_BarBIQ3, '--out', $outputfile_BarBIQ4, '--log', $log_file;
##Finished BarBIQ_sub_shift.pl##

##Run BarBIQ_sub_link.pl##
my $outputfile_BarBIQ42="$outputfile_BarBIQ4"."_link";
unlink $outputfile_BarBIQ42;
$code = $pippath.'BarBIQ_sub_link.pl';
     system $code, '--in',$outputfile_BarBIQ4,'--out',$outputfile_BarBIQ42, '--log', $log_file;
##Finished BarBIQ_sub_link.pl##

##Run BarBIQ_sub_link_error.pl##
my $outputfile_BarBIQ43="$outputfile_BarBIQ42"."_derr";
unlink $outputfile_BarBIQ43;
$code = $pippath.'BarBIQ_sub_link_error.pl';
system $code, '--in',$outputfile_BarBIQ42,'--out',$outputfile_BarBIQ43, '--log', $log_file;
#Finished BarBIQ_sub_link_error.pl##

##Run BarBIQ_sub_indels.pl##
my $outputfile_BarBIQ5="$outputfile_BarBIQ43"."_indels";
unlink $outputfile_BarBIQ5;
$code = $pippath.'BarBIQ_sub_indels.pl';
system $code, '--in',$outputfile_BarBIQ43,'--out',$outputfile_BarBIQ5, '--log', $log_file;
##Finished BarBIQ_sub_indels.pl##

##Run BarBIQ_sub_chimeras.pl##
my $outputfile_BarBIQ6="$outputfile_BarBIQ5"."_chimeras";
unlink $outputfile_BarBIQ6;
$code = $pippath.'BarBIQ_sub_chimeras.pl';
system $code, '--in',$outputfile_BarBIQ5,'--out',$outputfile_BarBIQ6, '--log', $log_file;
##Finished BarBIQ_sub_chimeras.pl##

##Run BarBIQ_sub_similar.pl##
my $outputfile_BarBIQ7="$outputfile_BarBIQ6"."_similar";
unlink $outputfile_BarBIQ7;
$code = $pippath.'BarBIQ_sub_similar.pl';
system $code, '--in',$outputfile_BarBIQ6,'--out',$outputfile_BarBIQ7, '--log', $log_file;
##Finished BarBIQ_sub_similar.pl##

##Run BarBIQ_sub_review_statistic.pl##
my $outputfile_BarBIQ8 = "$outputfile_BarBIQ7"."_statistic";
unlink $outputfile_BarBIQ8;
$code = $pippath.'BarBIQ_sub_statistic.pl';
system $code, '--in',$outputfile_BarBIQ7,'--out',$outputfile_BarBIQ8, '--log', $log_file;
##Finished BarBIQ_sub_statistic.pl##

##Run BarBIQ_sub_clustering_step_two.pl 0.1##
my $outputfile_BarBIQ01 = "$outputfile_BarBIQ2"."_0.1_sub2"; 
unlink $outputfile_BarBIQ01;
$code = $pippath.'BarBIQ_sub_clustering_step_two.pl';
system $code,'--in',$outputfile_BarBIQ2,'--out', $outputfile_BarBIQ01, "--ratio", "0.1", "--log", $log_file;
##Finished BarBIQ_sub_clustering_step_two.pl 0.1##

##Run BarBIQ_sub_shift.pl##
my $outputfile_BarBIQ02="$outputfile_BarBIQ01"."_shift";
unlink $outputfile_BarBIQ02;
$code = $pippath.'BarBIQ_sub_shift.pl';
system $code, '--in',$outputfile_BarBIQ01,'--out',$outputfile_BarBIQ02, '--log', $log_file;
##Finished BarBIQ_sub_shift.pl##

##Run BarBIQ_sub_link.pl##
my $outputfile_BarBIQ03="$outputfile_BarBIQ02"."_link";
unlink $outputfile_BarBIQ03;
$code = $pippath.'BarBIQ_sub_link.pl';
system $code, '--in',$outputfile_BarBIQ02,'--out',$outputfile_BarBIQ03, '--log', $log_file;
##Finished BarBIQ_sub_link.pl##

## delete the middle files ##
if($keep_middle eq "No")
     {
      print "You are going to delete the intermediate files\n";
      unlink ($outputfile_BarBIQ1, $outputfile_BarBIQ2, $outputfile_BarBIQ3, $outputfile_BarBIQ4, $outputfile_BarBIQ42, $outputfile_BarBIQ43, $outputfile_BarBIQ5, $outputfile_BarBIQ6);
      unlink ($outputfile_BarBIQ01, $outputfile_BarBIQ02);
      my $del_file1 = $outputfile_BarBIQ03."_bases";
      unlink $del_file1;
      my $del_file2 = $outputfile_BarBIQ42."_bases";
      unlink $del_file2;
      my $del_file3 = $outputfile_BarBIQ02."_statistic_I1";
      my $del_file4 = $outputfile_BarBIQ02."_statistic_R2";
      unlink ($del_file3, $del_file4);
      my $del_file5 = $outputfile_BarBIQ4."_statistic_I1";
      my $del_file6 = $outputfile_BarBIQ4."_statistic_R2";
      unlink ($del_file5, $del_file6);
      my $del_file7 = $outputfile_BarBIQ43."_errors";
      unlink $del_file7;
      my $del_file8 = $outputfile_BarBIQ01."_sort";
      unlink $del_file8;
      my $del_file9 = $outputfile_BarBIQ6."_cases_no";
      unlink $del_file9;
      my $del_file10 = $outputfile_BarBIQ7."_ratio";
      unlink $del_file10;
      my $del_file11 = $outputfile_BarBIQ5."_statistic";
      unlink $del_file11;
      my $del_file12 = $outputfile_BarBIQ3."_sort";
      unlink $del_file12;
      if ($inputfile_I1 =~ m/_I1_/s)
              {
                 unlink "$`_R1_$'";
                 unlink "$`_I1_$'";
                 unlink "$`_R2_$'";
              }
          else{die "Your input is wrong 0007777!!! please check!!! $!";}
     }    
## delete the middle files ##
print "$0 has done at:";
print scalar localtime;
print "\n";

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.003
#2023.02.16
