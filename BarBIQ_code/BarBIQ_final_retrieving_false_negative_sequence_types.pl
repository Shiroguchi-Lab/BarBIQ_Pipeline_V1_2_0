#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used to retrieve correct Bar sequences.
##########################################################################################################################################
#####how to run this code #####
##command##
#BarBIQ_final_review_overlap.pl --bar Bar_sequnence_file --in inputfile
##explaination##
#in: the input file which is the output file from BarBIQ_sub_link.pl and used parameter 0.1 in BarBIQ_sub_clustering_step_two.pl.
#Bar_sequnence_file: the output file from BarBIQ_final_library.pl.
##########################################################################################################################################
#####Install#####
## NONE
##########################################################################################################################################

#####code#####
use strict;
use warnings;

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";
##read command##
my ($i,$number,$inputfile,$outputfile, $stat_file, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--bar") {$stat_file = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --bar file --log log_file\"\n $!";}
    }
if(!$log_file) {
                 $log_file = $inputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
$outputfile="$inputfile"."_Rev";
unlink $outputfile;
if(!$stat_file)  {die "Your input is wrong!!!\n Please input \"--compare file\"\n $!";}
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
if(!(-e $stat_file)) {die "Your input file $stat_file is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfiles is:\n$inputfile\n";
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
if(!(($#info == 5) && ($info[0] =~ /\Acluster_/)))
   {
    die "Your input file $inputfile is wrong!!! please check!!!\n $!";
   }
close FILE;
open (FILE,$stat_file) or die "Could not open file '$stat_file' $!"; 
print "Your Bar-sequence data is:\n$stat_file\n";
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
if(!($info[0] eq "Bar-sequence-ID"))
   {
    die "Your input file $stat_file is wrong!!! please check!!!\n $!";
   }
close FILE;
print "Inputfiles are OK!\nStart to calculating:\n";
##check the inputfile##

###Main code### 
    print "Reviewing the original sequences...\n";
    my %keep_seq;
    open(STAT, $stat_file) or die "canot open input file '$stat_file' $!";
    $gi=<STAT>;
    @info=split(/\s+/,$gi);
    my $seq;
    for (my $i=0; $i<=$#info; $i++)
       {
        if ($info[$i] eq "Sequence") {$seq = $i;}
       }
    if (defined $seq) {} else {die "Your input file $stat_file is wrong!!! please check!!!\n $!";}
  #    my $seq = $#info-1;
        while($gi=<STAT>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);  
          $keep_seq{$info[$seq]}=$info[0]; 
         }
    close STAT;
    
    open (FILE,$inputfile) or die "Could not open file '$inputfile' $!";
    open(OUTF, '>>', $outputfile) or die "canot open output file '$outputfile' $!";
    while($gi=<FILE>)
          {
            chomp $gi;
            my @info=split(/\s+/, $gi);
            if(exists $keep_seq{$info[5]})
               {
                 print OUTF ("$gi\n");
               }
           }
    close FILE;
    close OUTF;
    print "Finished!!!\n";
print "Done!!!\n";
###Main code###
####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.19
