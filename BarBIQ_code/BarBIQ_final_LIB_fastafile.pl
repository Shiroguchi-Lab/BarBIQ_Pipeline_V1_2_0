#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code it used to generate a fasta file of the Bar sequences 
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_LIB_fastafile.pl --lib file1
##explaination##
#file1: the output file from BarBIQ_final_lib_COTU_clean.pl
####################################################################################################

#####code#######
use strict;
use warnings;

##save the input file names and output file name##
# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";
##read command##
my ($i,$libary,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--lib")  {$libary = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];} 
     else                         {die "Your input is wrong!!!\n Please input \"--lib Libary (--log log_file)\"\n $!";}
    }
if(!$libary)  {die "Your input is wrong!!!\n Please input \"--libary file\"\n $!";}

##read command##
##check the output name##
my $outputfile=$libary.".fa";
unlink $outputfile;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

     open(FILE, $libary) or die "cannot open input file '$libary' $!";
     open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
     my $gi=<FILE>;
     chomp $gi;
     my @info=split(/\s+/,$gi);
     my $seq;
     my $seqid;
     for(my $i=0; $i<=$#info; $i++)
         {
          if($info[$i] eq "Sequence") {$seq = $i;}
          if($info[$i] eq "Bar-sequence-ID") {$seqid = $i;}
         }
     if(defined $seq && defined $seqid) {} else{die "Your input is wrong!!!\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          print OUTF (">$info[$seqid]\n$info[$seq]\n");
         }
      close FILE; 
      close OUTF;
     
print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.001
##2018.12.14
