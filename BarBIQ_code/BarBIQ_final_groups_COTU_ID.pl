#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
##This code is used to add the cOTU IDs into the output file from BarBIQ_final_overlap_groups.pl
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_groups_COTU_ID.pl --group file1 --lib file2
##explaination##
#file1: output file from BarBIQ_final_overlap_groups.pl
#file2: the output file from BarBIQ_final_lib_COTU_ID.pl
########################################################################################################## 
#####Install#####
##None
############################################################################################################

#####code#######
use strict;
use warnings;

##save the input file names and output file name##
# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";
##read command##
my ($i,$groupfile, $libary,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--group")  {$groupfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--bar") {$libary = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}   
     else                         {die "Your input is wrong!!!\n Please input \"--group file --bar bar-sequence-file --log log_file\"\n $!";}
    }
if(!$groupfile)   {die "Your input is wrong!!!\n Please input \"--group file\"\n $!";}
if(!$libary)  {die "Your input is wrong!!!\n Please input \"--bar bar-sequence-file\"\n $!";}
##read command##
##check the output name##
my $outputfile=$groupfile."_COTU";
unlink $outputfile;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

my %COTU;
# my %Barseq;
     open(FILE, $libary) or die "cannot open input file '$libary' $!";
     my $gi=<FILE>;
     chomp $gi;
     my @info=split(/\s+/,$gi);
     my $COTU;
     my $seqid;
     my $orgid;
     for(my $i=0; $i<=$#info; $i++)
         {
          if($info[$i] eq "cOTU-ID") {$COTU = $i;}
          if($info[$i] eq "Bar-sequence-ID") {$seqid = $i;}
          if($info[$i] eq "OriginalID") {$orgid = $i;}
         }
     if(defined $COTU && defined $seqid && defined $orgid) {} else{die "Your input is wrong!!!\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          $COTU{$info[$seqid]}=$info[$COTU];
          # $Barseq{$info[$orgid]}=$info[$seqid];
         }
      close FILE; 
     
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
open(FILE1, $groupfile) or die "cannot open input file '$groupfile' $!";
while($gi=<FILE1>)
           {
             chomp $gi;
             @info=split(/\s+/,$gi);
             $info[0] = $COTU{$info[1]};
             # $info[1] = $Barseq{$info[1]};
             # $info[2] = $Barseq{$info[2]};
             my $out=join("\t", @info);
             print OUTF ("$out\n");
           }
close OUTF; 
close FILE1;

print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.20
