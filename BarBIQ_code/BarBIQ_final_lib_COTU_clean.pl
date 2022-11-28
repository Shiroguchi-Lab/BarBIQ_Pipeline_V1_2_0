#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to remove the contaminated cOTUs from the Bar sequences file
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_lib_COTU_clean.pl --COTULIB file1 --bcc file2
##explaination
#file1: the output file from BarBIQ_final_lib_COTU_ID.pl
#file2: the output file from BarBIQ_final_compare_datasets.pl
####################################################################################################

#####code#######
use strict;
use warnings;

##save the input file names and output file name##
# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";
##read command##
my ($i,$bac_count, $libary_COTU, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--COTULIB")  {$libary_COTU = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--bcc") {$bac_count = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--COTULIB libary_COTU --bcc bac_count --log log_file\"\n $!";}
    }

if(!$libary_COTU)   {die "Your input is wrong!!!\n Please input \"--COTULIB libary_COTU\"\n $!";}
if(!$bac_count)  {die "Your input is wrong!!!\n Please input \"--bcc bac_count\"\n $!";}
##read command##
##check the output name##
my $outputfile=$libary_COTU."_clean";
unlink $outputfile;
if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
my %COTU_kept;
     open(FILE, $bac_count) or die "cannot open input file '$bac_count' $!";
     my $gi=<FILE>;
     chomp $gi;
     my @info=split(/\s+/,$gi);
     if($info[0] eq "ID") { } else {die "Your input is wrong002!!!\n";}
     while($gi=<FILE>)
         {
           chomp $gi;
           @info=split(/\s+/,$gi);
           $COTU_kept{$info[0]} = $info[0];
         }
     close FILE;

     open(FILE, $libary_COTU) or die "cannot open input file '$libary_COTU' $!";
     open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
     $gi=<FILE>;
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $COTU;
     my $seqid;
     for(my $i=0; $i<=$#info; $i++)
         {
          if($info[$i] eq "cOTU-ID") {$COTU = $i;}
          if($info[$i] eq "Bar-sequence-ID") {$seqid = $i;}
         }
     if(defined $COTU && defined $seqid) {} else{die "Your input is wrong005!!!\n";}
     print OUTF ("$gi\n");
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          if(exists $COTU_kept{$info[$COTU]}) 
                 { print OUTF ("$gi\n"); }
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
