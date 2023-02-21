#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used to select the reads I1 and R2 according to the cleaned R1 file.
## The cleaned R1 file is generated from the code BarBIQ_sub_index_delete_leakage.pl.
## Please prepare the I1 reads file and R2 reads file, clustered R1 file for the same sample in advance.
## The input I1 and R2 files should be fastq format and with .fastq or .fastq.gz or .fq or .fq.gz
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_clean_I1R2_by_R1.pl --I1 I1.fastq --R2 R2.fastq --cluster R1file_cluster
###interpretation##
## --I1: I1 reads file, which is directly from illumina Miseq output, should be fastq format, end with .fastq or .fastq.gz or .fq or .fq.gz
## --R2: R2 reads file, which is directly from illumina Miseq output, should be fastq format, end with .fastq or .fastq.gz or .fq or .fq.gz. The order of all reads in I1 and R2 files should be the same, the files from miseq directly will be suitable
## --cluster: Clustered R1 reads file, which is the outputfile from BarBIQ_sub_index_delete_leakage.pl, please check this R1 file is from the same run of reads I1&R2
##########################################################################################################################################
#####Install#####
## please install the perl Module IPC::System::Simple qw(system) before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system);

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile_I1,$inputfile_R2,$outputfile_I1,$outputfile_R2, $cluster_R1, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--I1") {$inputfile_I1 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--R2") {$inputfile_R2 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--cluster") {$cluster_R1 = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else  {die "Your input is wrong!!!\n Please input \"--I1 I1.fastq --R2 R2.fastq --cluster R1file_cluster (--log log_file_name)\"\n $!";}
    }
if(!$inputfile_I1)   {die "Your input is wrong!!!\n Please input \"--I1: I1.fastq\"\n $!";}
if(!$inputfile_R2)   {die "Your input is wrong!!!\n Please input \"--R2: R2.fastq\"\n $!";}
if(!$cluster_R1)   {die "Your input is wrong!!!\n Please input \"--cluster R1file_cluster\"\n $!";}
if(!$log_file) {
                 $log_file = $cluster_R1."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
##read command##

##check the inputfile##
   if($inputfile_I1 =~ /.fastq.gz\z/)
       {
        print "$inputfile_I1 is unzipping...\n";
        my $inputfile_I1_unzip=$`.".fastq";
        system "gunzip -c $inputfile_I1 > $inputfile_I1_unzip"; print "unzip finished!!!\n";
        $inputfile_I1=$`.".fastq";
       }
elsif($inputfile_I1 =~ /.fq.gz\z/)
       {
        print "$inputfile_I1 is unzipping...\n";
        my $inputfile_I1_unzip=$`.".fq";
        system "gunzip -c $inputfile_I1 > $inputfile_I1_unzip"; print "unzipping finished!!!\n";
        $inputfile_I1=$`.".fq";
       }
elsif($inputfile_I1 =~ /.fastq\z/)
       {
       
       }
elsif($inputfile_I1 =~ /.fq\z/)
       {
        
       }
else   {die "Your inputfile '$inputfile_I1' is not fastq/fastq.gz/fq/fq.gz!!! please check!!! $!";}

if($inputfile_R2 =~ /.fastq.gz\z/)
       {
        print "$inputfile_R2 is unzipping...\n";
        my $inputfile_R2_unzip=$`.".fastq";
        system "gunzip -c $inputfile_R2 > $inputfile_R2_unzip"; print "unzip finished!!!\n";
        $inputfile_R2=$`.".fastq";
       }
elsif($inputfile_R2 =~ /.fq.gz\z/)
       {
        print "$inputfile_R2 is unzipping...\n";
        my $inputfile_R2_unzip=$`.".fq";
        system "gunzip -c $inputfile_R2 > $inputfile_R2_unzip"; print "unzipping finished!!!\n";
        $inputfile_R2=$`.".fq";
       }
elsif($inputfile_R2 =~ /.fastq\z/)
       {

       }
elsif($inputfile_R2 =~ /.fq\z/)
       {

       }
else   {die "Your inputfile '$inputfile_R2' is not fastq/fastq.gz/fq/fq.gz!!! please check!!! $!";}


my ($outputfile_I1_ave_qual, $outputfile_R2_ave_qual);
if($inputfile_I1 =~ m{.fastq\Z})
    {$outputfile_I1 = "$`"."_R1.fastq"; $outputfile_I1_ave_qual="$`"."_R1_ave_qual";}
else
    {die "Your inputfile '$inputfile_I1' is not fastq!!! please check!!! $!";}

if($inputfile_R2 =~ m{.fastq\Z})
    {$outputfile_R2 = "$`"."_R1.fastq"; $outputfile_R2_ave_qual="$`"."_R1_ave_qual";}
else
    {die "Your inputfile '$inputfile_R2' is not fastq!!! please check!!! $!";}
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
##check the inputfile##

##find the total reads for printing the progress##
open (FILE,$inputfile_I1) or die "Could not open file '$inputfile_I1' $!"; # open inputfile again, read after the last line.
$gi=<FILE>;$gi=<FILE>;$gi=<FILE>;$gi=<FILE>;
my $position_1=tell(FILE);
seek (FILE, 0, 2);
my $position_2=tell(FILE);
my $reads=$position_2/$position_1/10;
close FILE;
##find the total reads for printing the progress##

##read the cluster informations##
open (CLUS,$cluster_R1) or die "Could not open file '$cluster_R1' $!";
my %clusters;
while($gi=<CLUS>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    $clusters{$info[0]} = $info[3];
   }
close CLUS;
##read the cluster informations##

##Main process##
open (FILEI1,$inputfile_I1) or die "Could not open file '$inputfile_I1' $!";
open (FILER2,$inputfile_R2) or die "Could not open file '$inputfile_R2' $!";
open (OUTFI1,'>>', $outputfile_I1) or die "Could not open file '$outputfile_I1' $!";
open (OUTFR2,'>>', $outputfile_R2) or die "Could not open file '$outputfile_R2' $!";
my ($gi_I1, $gi_R2,@info_I1,@info_R2,$seq_I1,$seq_R2,$qual_I1,$qual_R2,$ID_I1,$ID_R2);
my $No_reads=0;
my $unbarcode_no=0;
my $p=10;
my $progress=$reads; 
while($gi_I1=<FILEI1>)
   {
     $gi_R2=<FILER2>;
     $No_reads++; if($No_reads>=$progress) { print "$p"; print"%\n"; $progress=$progress+$reads;$p=$p+10;} # print the progress
     chomp $gi_I1;
     chomp $gi_R2;
     @info_I1=split(/\s+/,$gi_I1); 
     @info_R2=split(/\s+/,$gi_R2);
     if($info_I1[0] =~ m/\A@(.*)/s) # read the illumina ID of I1
         {
          $ID_I1=$1;
          if($info_R2[0] =~ m/\A@(.*)/s) # read the illumina ID of R2
              {
               $ID_R2=$1;
               if ($ID_I1 ne $ID_R2) {die "Your I1&R2 files have different orders of reads\n $!";}
              }
          else{die "Your inputfile '$inputfile_R2' is not fastq format, please check\n $!"}
         }
     else{die "Your inputfile '$inputfile_I1' is not fastq format, please check\n $!"}
     chomp($gi_I1=<FILEI1>); $seq_I1=$gi_I1;# read the sequence of I1
     chomp($gi_R2=<FILER2>); $seq_R2=$gi_R2;# read the sequence of R2
     my $I1=<FILEI1>; my $R2=<FILER2>;
     chomp($gi_I1=<FILEI1>); $qual_I1=$gi_I1;# read the quality of I1
     chomp($gi_R2=<FILER2>); $qual_R2=$gi_R2;# read the quality of R2
     if(exists $clusters{$ID_I1})
        {
         print OUTFI1 ("\@$ID_I1\n$seq_I1\n$I1$qual_I1\n");
         print OUTFR2 ("\@$ID_R2\n$seq_R2\n$R2$qual_R2\n");
         }
      else {$unbarcode_no++}
     }
close FILEI1;
close FILER2;
close OUTFI1;
close OUTFR2; 
my $barcode_no=$No_reads-$unbarcode_no;
my $kep_ratio=$barcode_no/$No_reads; # The ratio of kept reads
print "The No. of unbarcoded reads: $unbarcode_no\nThe No. of barcoded reads: $barcode_no\nThe ratio of barcoded reads: $kep_ratio\n";

print "Done!!!\n";
##Main process##

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
