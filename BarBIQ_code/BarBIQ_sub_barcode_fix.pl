#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used for deleting the barcode according to the fixed base which were designed in the original barcode
## You should prepare a file which contain the designed fixed bases like the BarBIQ_example_fixed_base.txt and for the real sample and standard samples individuely
## You can do at the same time for mutipule files which used the same fixed bases. Finally, you will get a merged file for all files
## Please prepare a file which contain the file names of all inputfiles like the BarBIQ_example_inputfile_name.txt
## The input data should be fasta format and with .fasta, output file will also be fasta format.
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_barcode_fix.pl --in BarBIQ_example_inputfile_name.txt --out outputfile --fixed BarBIQ_example_fixed_base.txt
###interpretation##
## --in: a file which contain the file names of all inputfiles, should be prepared like BarBIQ_example_inputfile_name.txt 
## --out: outputfile, please set a name for your outputfile, should end with .fasta
## --fixed: a file which contain the designed fixed bases, should be prepared like the BarBIQ_example_fixed_base.txt
##########################################################################################################################################
#####Install#####
## Please install the perl model: Bio::SeqIO and Bio::Seq before use this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use Bio::SeqIO; # for reading the fasta files please install them before run this code
use Bio::Seq;

# print "Now you are runing BarBIQ_sub_barcode_fix.pl\n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$outputfile,$fixed_base, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--fixed") {$fixed_base = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];} 
     else                         {die "Your input is wrong!!!\n Please input \"--in example_inputfile_name.txt --out outputfile --fixed example_fixed_base.txt (--log log_file_name)\"\n $!";}
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

if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
if(!$fixed_base)  {die "Your input is wrong!!!\n Please input \"--fixed example_fixed_base.txt\"\n $!";}
##read command##

##check the inputfile##
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
my (@inputfiles, $gi, @Fixed_bases);
print "Your inputfiles are:\n";
while($gi=<FILE>)
    {
     chomp $gi;
     if($gi =~ /.fasta/)
         {
           my $file_name="$`".".fasta";
           push @inputfiles, $file_name;
           print "$file_name\n";
         }
    }
close FILE;
open (FIXED, $fixed_base) or die "Could not open file '$fixed_base' $!"; # open fixed_base file
print "Your fixed bases are:\n";
while($gi=<FIXED>)
    {
     chomp $gi;print "$gi\n";
     push @Fixed_bases, $gi;
    }
close FIXED;

for($i=0; $i<=$#inputfiles; $i++)
    {
     if(!(-e $inputfiles[$i])) {die "Your input file '$inputfiles[$i]' is not existed!!! please check!!!\n $!";}
     else {open (FILEFASTA,$inputfiles[$i]) or die "Could not open file '$inputfiles[$i]' $!";
           chomp ($gi=<FILEFASTA>);
           if($gi=~m/\A>(.*)/s)
               {
                $gi=<FILEFASTA>; chomp ($gi=<FILEFASTA>);
                if(!($gi=~m/\A>(.*)/s)) {die "Your input file $inputfiles[$i] is not fasta format!!! please check!!!\n $!";}
               }
           else {die "Your input file $inputfiles[$i] is not fasta format!!! please check!!!\n $!";}
           close FILEFASTA;
          }
    }
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

##find the total reads for printing the progress##
my $reads=0;
for($i=0; $i<=$#inputfiles; $i++)
    {
     open (FILE,$inputfiles[$i]) or die "Could not open file '$inputfiles[$i]' $!";
     $gi=<FILE>;$gi=<FILE>;
     my $position_1=tell(FILE);
     seek (FILE, 0, 2);
     my $position_2=tell(FILE);
     $reads=$reads+$position_2/$position_1;
    }
$reads=$reads/10;
close FILE;
##find the total reads for printing the progress##

##check each read and delete the unmatched read to the fixed bases##
my ($display_name, $desc, $DNA);
my $No=0; # The No. of total reads
my $keep_read=0; # The No. of ketp reads
my $progress=$reads; # for printing the progress
my $p=10; # for printing the progress
open(OUTF, '>>', $outputfile) or die "Could not open file '$outputfile' $!";
for($i=0; $i<=$#inputfiles; $i++)
     {
      my $catchseq_seqio_obj = Bio::SeqIO -> new(-file   => $inputfiles[$i], #read the fasta file
                                                 -format => 'fasta');
      while(my $seq_obj = $catchseq_seqio_obj->next_seq)
          {
           $No++;
           if($No>=$progress) { print "$p"; print"%\n"; $progress=$progress+$reads;$p=$p+10;} # print the progress 
           $display_name = $seq_obj->display_name; # read the display name of each read
           $desc = $seq_obj->desc;  # read the description of each read
           $DNA= $seq_obj->seq; # read the DNA sequence of each read
           &Fixed_base_check($display_name, $desc, $DNA);
          }
     }
close OUTF;

my $del_read=$No-$keep_read; # The No. of deleted reads
my $kep_ratio=$keep_read/$No; # The ratio of kept reads
print "The No. of deleted reads: $del_read\nThe No. of kept reads: $keep_read\nThe ratio of kept reads: $kep_ratio\n";
print "Done!!!\n";

sub Fixed_base_check # check the read's match to the fixed bases or not, and decide to keep or delete
    {
     my ($j, $delete);
     $delete=0;
     for($j=0; $j<=$#Fixed_bases; $j++)
         {
           if($_[2]=~ m{$Fixed_bases[$j]}) {$delete=1; last;}
         }
     if($delete==1)
         {
          print OUTF (">$_[0] $_[1]\n$_[2]\n");
          $keep_read++;
         }
    }
##check each read and delete the unmatched read to the fixed bases##

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.003
#2022.10.14
