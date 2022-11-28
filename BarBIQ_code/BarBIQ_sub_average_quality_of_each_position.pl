#! /usr/bin/env perl
##########################################################################################################
#####Description of this code#####
# This code is used to determine the trimming length for each sequencing run.
# Average quality scores for two consecutive bases will be calculated for all positions and all reads. 
# The trimming length is use to trim the 3’ ends of both I1 and R2 reads, which tend to contain errors due to their low sequencing qualities.
##########################################################################################################
#####how to run this code #####
##command##
#BarBIQ_sub_average_quality_of_each_position.pl --in file --out file
##interpretation##
# --in:  inputfile, which has to be fastq file, file.fastq file.fq or file.fastq.gz file.fq.gz are accepted.
# --cut: filename, please ste a name for your outputfile to save the suggested cut postion
# --out: outputfile, please set a name for your outputfile
###########################################################################################################
#####Install#####
## please install the perl Module IPC::System::Simple qw(system) before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system); # for calling external program

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$cut,$outputfile, $log_file); 
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--cut") {$cut = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --cut filename --out: outputfile (--log log_file_name)\" $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in inputfile\" $!";}
if(!$cut)   {die "Your input is wrong!!!\n Please input \"--cut filename\" $!";}
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out outputfile\" $!";}
if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!! $!";}
##read command##

##check the inputfile##
my ($inputfile_fastq,$gi,@info);
   if($inputfile =~ /.fastq.gz\z/) 
       {
        print "$inputfile is unzipping...\n"; 
        $inputfile_fastq=$`.".fastq";
        system "gunzip -c $inputfile > $inputfile_fastq"; print "unzip finished!!!\n"; 
       }
elsif($inputfile =~ /.fq.gz\z/)
       {
        print "$inputfile is unzipping...\n";
        $inputfile_fastq=$`.".fq";
        system "gunzip -c $inputfile > $inputfile_fastq"; print "unzipping finished!!!\n";
       }
elsif($inputfile =~ /.fastq\z/)
       {
        $inputfile_fastq=$inputfile;
       }
elsif($inputfile =~ /.fq\z/)
       {
        $inputfile_fastq=$inputfile;
       }
else   {die "Your inputfile '$inputfile' is not fastq/fastq.gz/fq/fq.gz!!! please check!!! $!";}
open (FILE,$inputfile_fastq) or die "Could not open file '$inputfile_fastq' $!"; # open inputfile
chomp ($gi=<FILE>);
   if($gi =~ m/\A@(.*)/s) # check the first line
       {
        $gi=<FILE>;$gi=<FILE>;$gi=<FILE>;chomp ($gi=<FILE>); # read to the 5th line
        if($gi =~ m/\A@(.*)/s) # check the 5th line
            {print "Inputfile is OK!\nStart to calculating:\n";}
      else  {die "Your inputfile '$inputfile' is not fastq!!! please check!!! $!";}
       }
 else  {die "Your inputfile '$inputfile' is not fastq format!!! please check!!! $!";}
close FILE; # close inputfile
##check the inputfile##

##find the total reads for printing the progress##
open (FILE,$inputfile_fastq) or die "Could not open file '$inputfile_fastq' $!"; # open inputfile again, read after the last line.
$gi=<FILE>;$gi=<FILE>;$gi=<FILE>;$gi=<FILE>;
my $position_1=tell(FILE);
seek (FILE, 0, 2);
my $position_2=tell(FILE);
my $reads=$position_2/$position_1/10;
close FILE;
##find the total readsfor printing the progress##

##calculate the average quality##
open (FILE,$inputfile_fastq) or die "Could not open file '$inputfile_fastq' $!"; # open inputfile again, read from the first line.
my ($No, $line, $quality, @quality, @sum_quality, @ave_quality, $average_qual, $pos, $recommend_cut);
my $p=10;
my $progress=$reads;
chomp ($gi=<FILE>);
   if ($gi =~ m/\A@(.*)/s)
       {
        $gi=<FILE>;
        $gi=<FILE>;
        chomp ($gi=<FILE>); $quality = $gi; # read the first read's quality  
        @quality=split(//,$quality);
        for ($i=0; $i<=$#quality; $i++)
            {
             $sum_quality[$i] = ord($quality[$i])-33; #calculate the sum of qualities for each position
            }
        $No=1;
       }
 else  {$line=($No-1)*4+1; die "Your inputfile '$inputfile' is not fastq format at line $line!!! please check!!! $!";}
while($gi=<FILE>)
     {
      chomp $gi;
      if ($gi =~ m/\A@(.*)/s)
          {
           $gi=<FILE>;
           $gi=<FILE>;
           chomp ($gi=<FILE>); $quality = $gi; # read quality
           @quality=split(//,$quality);
           for ($i=0; $i<=$#quality; $i++)
               {
                $sum_quality[$i] = $sum_quality[$i]+ord($quality[$i])-33; #calculate the sum of qualities for each position
               }
           $No++; if($No>=$progress) { print "$p"; print"%\n"; $progress=$progress+$reads;$p=$p+10;} # print the progress
          }
      else  {$line=($No-1)*4+1; die "Your inputfile '$inputfile' is not fastq format at line $line!!! please check!!! $!";}
     }
open (OUTF,'>>', $outputfile) or die "Could not open file '$outputfile' $!"; # open outputfile as a writing file
for ($i=0; $i<=$#sum_quality; $i++)
    {
     $ave_quality[$i] = $sum_quality[$i]/$No; #calculate the average qualities for each position
     $pos=$i+1;
     print OUTF ("$pos\t$ave_quality[$i]\n"); #print the average qualities for each position
    }
$recommend_cut = $#ave_quality+1; ## if all positions are more than 25, keep all.
for ($i=5; $i<=$#ave_quality; $i++)
    {
     my $ave_quality_ave = ($ave_quality[$i]+$ave_quality[$i-1])/2;
     if ($ave_quality_ave < 25) { $recommend_cut = $i; last;} #The recommended cut off position
    }
close FILE;
close OUTF;
my $outputfile2 = $cut;
unlink $outputfile2;
open (OUTF,'>>', $outputfile2) or die "Could not open file2 $!"; # open outputfile as a writing file
print OUTF ("$recommend_cut\n");
close OUTF;
print "Done!!!\n";
##calculate the average quality##

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.003
#2022.10.13
