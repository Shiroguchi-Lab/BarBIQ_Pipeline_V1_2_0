#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used for delete the leakage and seperate reads into individual files for each index.
## Leakage means different indexes (sample) have been clustered into a same cluster, this code considers the minority index as leakage and delete them.
## If the No. of reads for both indexes are the same, this code will delete both.
## After leakage deletion, all clustered reads will be seperated into individual files according to the indexes. 
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_index_and_leakage.pl --in inputfile --out outputfile --index BarBIQ_example_index.txt
###interpretation##
## --in: the outputfile from BarBIQ_sub_barcode_clustering, which is clustered by nucleotide-sequence-clusterizer
## --out: outputfile, please set a name for your outputfile, which will be used for generating series outputfiles
## --index: a file include all indexes for your inputfile, please prepare it like the given example, BarBIQ_example_index.txt
##########################################################################################################################################
#####Install#####
## None.
##########################################################################################################################################

#####code#####
use strict;
use warnings;

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$outputfile,$index,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--index") {$index = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in XXXX --out XXXXX --index XXXX (--log log_file_name)\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out: outputfile\"\n $!";}
if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

if(!$index)  {die "Your input is wrong!!!\n Please input \"--index BarBIQ_example_index.txt\"\n $!";}
##read command##

##check the inputfile##
my $gi;
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
$gi=<FILE>;
chomp $gi;
if ($gi =~ m{\A>([0-9]*)}) 
    {
      chomp ($gi=<FILE>);
      if ($gi =~ m{\A>([0-9]*)}) {die "Your input is wrong!!!\n The inputfile is not the outputfile of nucleotide-sequence-clusterizer\n $!";}
    }
else {die "Your input is wrong!!!\n The inputfile is not the outputfile of nucleotide-sequence-clusterizer\n $!";}
close FILE;
my @index;
open (FILE,$index) or die "Could not open file '$index' $!"; # open inputfile
chomp ($gi=<FILE>);
@index=split(/\s+/, $gi);
if ($index[0] ne "indexes"){die "Your input is wrong!!!\n The index file has not followed the BarBIQ_example_index.txt\n $!";}
print "Your indexes are: @index\n";
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

##find the total reads for printing the progress##
my $ttreads=0;
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!";
$gi=<FILE>;$gi=<FILE>;
my $position_1=tell(FILE);
seek (FILE, 0, 2);
my $position_2=tell(FILE);
$ttreads=$position_2/$position_1/10;
close FILE;
##find the total reads for printing the progress##

##Open files for each index##
my (@outputfile, %outf);
for ($i=1; $i<=$#index; $i++) {
   $outputfile[$i] = "$outputfile"."_$index[$i]";  ##outputfile of $i index
   if(-e $outputfile[$i]){die "Your output file $outputfile[$i] is already existed!!! please check!!!\n $!";}
   open(my $hd, '>>', $outputfile[$i]) or die "Could not open file '$outputfile[$i]' $!";
   $outf{$index[$i]} = $hd; ## file handles of $i index
   }
##Open files for each index##

##read clusters and delete the leakage and seperate all indexes##
open(FILE, $inputfile) or die "Could not open file '$inputfile' $!"; ## open the inputfile
my $progress=$ttreads; # for printing the progress
my $p=10; # for printing the progress
my ($cluster_name, @reads, @info, $index_keep, $No_reads, $No_leakage, $No_kept, $ratio_kept);
$No_reads=0;
$No_kept=0;
chomp ($gi=<FILE>);
if ($gi =~ m{\A>([0-9]*)}) 
    {
     $cluster_name = "cluster_"."$1";
    }
else {die "Your input is wrong!!!\n The inputfile is not the outputfile of nucleotide-sequence-clusterizer\n $!";}
while($gi=<FILE>)
    {
     chomp $gi;
     if ($gi =~ m{\A>([0-9]*)})
         {
          $index_keep = &del_leakage;
          &seperate_index($index_keep, $cluster_name);
          $cluster_name = "cluster_"."$1";
          undef @reads;
         }
     else{
          @info=split(/\s+/,$gi);
          push @reads, [@info];   
          $No_reads++; if($No_reads>=$progress) { print "$p"; print"%\n"; $progress=$progress+$ttreads;$p=$p+10;} # print the progress
         }
     }
          $index_keep = &del_leakage;
          &seperate_index($index_keep, $cluster_name);    
for ($i=1; $i<=$#index; $i++) 
    {
     close $outf{$index[$i]};
    }
$No_leakage=$No_reads-$No_kept; # The No. of deleted reads
$ratio_kept=$No_kept/$No_reads;  # The ratio of kept reads
print "The No. of leakage reads: $No_leakage\nThe No. of kept reads: $No_kept\nThe ratio of kept reads: $ratio_kept\n";

##delete the leakage##
sub del_leakage {
my @reads_sort = sort{$a->[1] cmp $b->[1]} @reads;
my ($j, @No, $No, $index, $index_large);
$index=$reads_sort[0][1];
$No=1;
if($#reads_sort>=1)
    {
     for($j=1; $j<=$#reads_sort; $j++)
         {
          if($reads_sort[$j][1] eq $index)
              {
                $No++;
              }
          else{
                push @No, [($index, $No)];
                $index=$reads_sort[$j][1];
                $No=1;
              }
         }
     push @No, [($index, $No)];
    }
else{push @No, [($index, $No)];}
my @No_sort= sort{$b->[1] <=> $a->[1]} @No;
if ($#No_sort == 0)
    {
     $index_large = $No_sort[0][0];
    }
else{
     if($No_sort[0][1] > $No_sort[1][1])
         {
           $index_large = $No_sort[0][0];
         }
     else{$index_large="*";}
    }
return $index_large;
}
##delete the leakage##

##seperate each index to individual file##
sub seperate_index {
my $k;
for($k=0; $k<=$#reads; $k++)
    {
     if($reads[$k][1] eq $_[0])
         {
          print {$outf{$_[0]}} ("$reads[$k][0]\t$reads[$k][1]\t$reads[$k][2]\t$_[1]\n");
          $No_kept++;
         }
    }
}
##seperate each index to individual file##
##read clusters and delete the leakage and seperate all indexes##

print "Done!!!\n";

##end##
#####code#####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
