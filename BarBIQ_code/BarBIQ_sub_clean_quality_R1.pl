#! /usr/bin/env perl
##########################################################################################################
#####Description of this code#####
# This code is used to remove low-quality reads for the Barcodes, in our case is the R1 reads.
# You can choose different thresholds by parameters, --AveWindow X and --QualThreshold Y, which used to remove the reads: at least one window of X continuous bases with an average sequencing quality score (determined by MiSeq) less than Y. Default is X=4 and Y=15.
# The input data should be fastq format, output file will be fasta format and will add the index ID of this data to each read.
# How many reads deleted using this step will be shown on your standard output.
##########################################################################################################
#####how to run this code #####
##command##
# BarBIQ_sub_clean_quality_R1.pl --in inputfile --out outputfile --index SX (--AveWindow X --QualThreshold Y)
##interpretation##
# --in:  inputfile, which has to be fastq file, file.fastq file.fq or file.fastq.gz file.fq.gz are accepted.
# --out: outputfile, please set a name for your outputfile, should end with .fasta
# --index: the index ID of this data
# --AveWindow: the window length for caclulate the average quality, you can change by command --AveWindow, default value is 4
# --QualThreshold: the average quality threshold for removing the read, you can change by command --QualThreshold, default value is 15
###########################################################################################################
#####Install#####
## Please install the perl modules, IPC::System::Simple qw(system), before use this code
############################################################################################################
#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system); # for calling external program

# print "Now you are runing BarBIQ_sub_clean_quality_R1.pl\n";
# print "The parameters are: @ARGV\n";
##read command##
my ($i,$inputfile,$outputfile,$index, $log_file); 
my $ave_seq=4; # the window length for caclulate the average quality, you can change by command --AveWindow, default value is 4
my $qual_thres=15; # the quality threshold for delete the read, which used to be compared with the average quality of each window length sequence of each read, you can change by command --QualThreshold, default value is 15
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--index") {$index = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--AveWindow") {$ave_seq = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--QualThreshold") {$qual_thres = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}     
     else                         {die "Your input is wrong!!!\n Please input \"--in: inputfile, --out: outputfile and --index SX (--log log_file_name)\"\n $!";}
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
print "You are using the threshold of $ave_seq as the window length for caclulate the average quality\n";
print "You are using the quality threshold of $qual_thres to delete the reads\n";

if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
if(!$index)  {die "Your input is wrong!!!\n Please input \"--index SX\"\n $!";}
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
##find the total reads for printing the progress##

##check each read and delete the low quality one##
open (FILE,$inputfile_fastq) or die "Could not open file '$inputfile_fastq' $!"; # open inputfile again, read from the first line.
open (OUTF,'>>', $outputfile) or die "Could not open file '$outputfile' $!"; # open outputfile as a writing file
my ($No, $line, @read_name, $read_name, $seq, $quality);
my $p=10;
my $progress=$reads;
$No=0; # The No. of total reads
my $keep_read=0; # The No. of kept reads
while($gi=<FILE>)
     {
      chomp $gi;
      if ($gi =~ m/\A@(.*)/s)
          {
           @read_name=split(/\s+/,$1);
           $read_name=$read_name[0];
           chomp ($gi=<FILE>); $seq=$gi; # read sequence
           $gi=<FILE>;
           chomp ($gi=<FILE>); $quality = $gi; # read quality
           &Low_quality_check($seq, $quality); #  check the read's quality and decide to keep or delete
           $No++; if($No>=$progress) { print "$p"; print"%\n"; $progress=$progress+$reads;$p=$p+10;} # print the progress
          }
      else  {$line=($No)*4+1; die "Your inputfile '$inputfile' is not fastq format at line $line!!! please check!!! $!";}
     }
close FILE;
close OUTF;
my $del_read=$No-$keep_read; # The No. of deleted reads
my $kep_ratio=$keep_read/$No; # The ratio of kept reads
print "The No. of deleted reads: $del_read\nThe No. of kept reads: $keep_read\nThe ratio of kept reads: $kep_ratio\n";
print "Done!!!\n";

sub Low_quality_check # check the read's quality and decide to keep or delete
    {
     my ($j, $k,  $sequence, @quality, $sum_quality, $average_quality, $delete);
     $sequence= $_[0];
     @quality=split(//,$_[1]);
     $delete=0;
     for($j=0; $j<=($#quality-$ave_seq+1); $j++)
         {
          $sum_quality=0;
          for($k=$j; $k<=($j+$ave_seq-1); $k++)
              {
               $sum_quality=$sum_quality+ord($quality[$k])-33;
              }
          $average_quality=$sum_quality/$ave_seq;
          if($average_quality<$qual_thres)
              {
               $delete=1;
               last;
              }
         }
     if($delete==0)
         {
          print OUTF (">$read_name $index\n$sequence\n");
          $keep_read++;
         }
    }
##check each read and delete the low quality one##

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2018.08.06
