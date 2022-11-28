#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to format the Bar sequence as a BarBIQ standard Bar sequence and cOTU file, which include the following information
#Bar-sequence-ID	ID for each Bar sequence.
#cOTU-ID	ID for each cellular-based operational taxonomy unit (cOTU).
#Seqeunce	Seqeunce of Bar sequence containning both I1 and R2 and linked by their overlopped sequences.
#Seqeunce-I1	Seqeunce of Bar sequence identified by I1 reads.
#Seqeunce-R2	Seqeunce of Bar sequence identified by R2 reads.
#LINK	If I1 and R2 are overlapped, labeled as "LINK" ; if not,  as "I1R2".
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_libary.pl inputfile log_file
##explaination##
#inputfile: the inputfile which is the output file from BarBIQ_final_merge_all_samples.pl
#log_file: a file to save the log information output by this code
####################################################################################################

#####code#######
use strict;
use warnings;

my $log_file = $ARGV[2];
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are running program: $0\n";
print "Your parameters are: @ARGV\n";

##save the input file names and output file name##
my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];
unlink $outputfile;
##check the output name##

##Main##
my %lib;
open(FILE, $inputfile) or die "cannot open input file '$inputfile' $!";
my $gi=<FILE>;
chomp $gi;my @info=split(/\s+/,$gi);
my ($id,$seq, $link, $I1, $R2);
my $no=1;
for(my $i=0; $i<=$#info; $i++)
    {
       if ($info[$i] eq "Sequence") {$seq = $i;}
    elsif ($info[$i] eq "LINK") {$link = $i;}
    elsif ($info[$i] eq "I1") {$I1 = $i;}
    elsif ($info[$i] eq "R2") {$R2 = $i;}
    elsif ($info[$i] eq "ID") {$id = $i;}
    }
if(defined $seq && defined $link && defined $I1 && defined $R2 && defined $id) {} else{die "Your inputfile is wrong!!\n"}
while($gi=<FILE>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    $lib{$no}=$gi;
    $no++;
   }
close FILE;

open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
print OUTF "Bar-sequence-ID\tcOTU-ID\tSequence\tI1\tR2\tLINK\tOriginalID\n";
my $bac_no="XXXXX";
foreach my $key (sort {$a <=> $b} keys %lib)
   {
     my @p=split(/\s+/,$lib{$key});
     my $pp=sprintf("%04s", $key);
     print OUTF ("Bar-sequence-$pp\t$bac_no\t$p[$seq]\t$p[$I1]\t$p[$R2]\t$p[$link]\t$p[$id]\n");
   }

close OUTF;
##end##

##main end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.18
