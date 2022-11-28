#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used for clusterizering the barcodes by a software nucleotide-sequence-clusterizer.
## nucleotide-sequence-clusterizer can be downloaded from http://kirill-kryukov.com/study/tools/nucleotide-sequence-clusterizer/.
## You should do at the same time for all indexes including samples and spike-in controls in a single MiSeq run. Finally, you will get clustered reads for each sample.
## Please prepare a file which contain the file names of all inputfiles like the example we provide, BarBIQ_example_inputfile_name.txt.
## The input data should be fasta format and with end of .fasta.
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_barcode_clustering.pl --in BarBIQ_example_inputfile_name.txt --out outputfile --d D --t T
###interpretation##
## --in: a file which contain the file names of all inputfiles, should be prepared like BarBIQ_example_inputfile_name.txt 
## --out: outputfile, please set a name for your outputfile
## --d: sequences are separated by D or fewer substitutions (a parameter for nucleotide-sequence-clusterizer, largest is 3, we suggest to use 2, however, 1, 2, and 3 do not change the results much)
## --t: No. of bases in the barcode, which is used as parameter -t for nucleotide-sequence-clusterizer
##########################################################################################################################################
#####Install#####
## Please install the nucleotide-sequence-clusterizer and add it to Environment variable (let it possible to be called directly)
## please install the perl Module IPC::System::Simple qw(system) before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system); ## for calling external program

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n"; 
##read command##
my ($i,$inputfile,$outputfile,$d,$t,$log_file, $software);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--soft") {$software = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--d") {$d = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--t") {$t = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];} 
     else                         {die "Your input is wrong!!!\n Please input \"--in example_inputfile_name.txt --out outputfile --soft nucleotide-sequence-clusterizer --d 2 --t 30 (--log log_file_name)\"\n $!";}
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
if(!$software) {die "Your input is wrong!!!\n Please input \"--soft nucleotide-sequence-clusterizer\"\n $!";}
if(!$d)  {die "Your input is wrong!!!\n Please input \"--d 2\"\n $!";}
if(!$t)  {die "Your input is wrong!!!\n Please input \"--t 30\"\n $!";}
##read command##

##check the inputfile##
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
my (@inputfile, $gi);
print "Your inputfiles are:\n";
while($gi=<FILE>)
    {
     chomp $gi;
     if($gi =~ /.fasta\Z/)
        {
          push @inputfile, $gi;
          print "$gi\n";
        }
    }
close FILE;

for($i=0; $i<=$#inputfile; $i++)
    {
     if(!(-e $inputfile[$i])) {die "Your input file $inputfile[$i] is not existed!!! please check!!!\n $!";}
     else {open (FILEFASTA,$inputfile[$i]) or die "Could not open file '$inputfile[$i]' $!";
           chomp ($gi=<FILEFASTA>);
           if($gi=~m/\A>(.*)/s)
               {
                $gi=<FILEFASTA>; chomp ($gi=<FILEFASTA>);
                if(!($gi=~m/\A>(.*)/s)) {die "Your input file $inputfile[$i] is not fasta format!!! please check!!!\n $!";}
               }
           else {die "Your input file $inputfile[$i] is not fasta format!!! please check!!!\n $!";}
           close FILEFASTA;
          }
    }
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

##merge all files##
print "Merging files...\n";
my $inputfile_merge="$outputfile"."_all_middle.fasta"; ## a file include all data of all inputfiles
unlink $inputfile_merge;
open (MERGE, '>>', $inputfile_merge) or die "Could not open file '$inputfile_merge' $!";
for($i=0; $i<=$#inputfile; $i++)
    {
     open (FILEFASTA,$inputfile[$i]) or die "Could not open file '$inputfile[$i]' $!";
     while($gi=<FILEFASTA>)
          {
           chomp $gi;
           print MERGE ("$gi\n");
          }
     close FILEFASTA;
    }
close MERGE;     
print "File Merging finished!!!\n";
##merge all files##

##nucleotide-sequence-clusterizer##
my $T="";
for($i=0; $i<$t; $i++)
   {$T="$T".".";}
system $software, '-i', $inputfile_merge, '-o', $outputfile, '-d', $d, '-t', $T; ## clustering using nucleotide-sequence-clusterizer
##nucleotide-sequence-clusterizer##
unlink $inputfile_merge;
print "Done!!!\n";

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
