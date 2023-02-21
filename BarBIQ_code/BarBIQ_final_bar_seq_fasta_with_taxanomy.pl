#! /usr/bin/env perl
###############################################################
#######Description of this code#####
#This code is used to write the Bar sequences into a fasta file with predicted taxanomy from RDP classifier.
##########################################################################
######how to run this code #####
##command##
#BarBIQ_final_bar_seq_fasta_with_taxanomy.pl --bar file1 --taxa file2 --out file3
##explaination##
#file1: the output file from  
#file2: the output file from BarBIQ_final_add_Taxonomy_to_bac_count_publish.pl 
#file3: outputfile.xlsx
########################################################################################################
#####Install#####
#Please install the module Excel::Writer::XLSX before using this code
##############################################################################################

#####code#####

use strict;
use warnings;
use Excel::Writer::XLSX;

##read command##
my ($i,$barseq_data,$outputfile,$taxanomy,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--bar") {$barseq_data = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--taxa") {$taxanomy = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
#     elsif ($ARGV[$i] eq "--dataname") {$dtnm = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--bar Bar_seq_COTU_clean --taxa All_samples_cell_count_annotation_RDP_classifier.txt --out outputfile (--log log_file)\"\n $!";}   
    }
if(!$barseq_data)  {die "Your input is wrong!!!\n Please input \"--bar Bar_seq_COTU_clean\"\n $!";}
# if(!$dtnm)  {die "Your input is wrong!!!\n Please input \"--dataname XX\"\n $!";}
# if(!$raw_data)  {die "Your input is wrong!!!\n Please input \"--raw raw_file\"\n $!";}
# if(!$index)  {die "Your input is wrong!!!\n Please input \"--index index_file\"\n $!";}
if(!$taxanomy)  {die "Your input is wrong!!!\n Please input \"--taxa All_samples_cell_count_annotation_RDP_classifier.txt\"\n $!";}
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out outputfile\"\n $!";}
if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
##read command##

##Taxanomy##
open(TAXA, $taxanomy) or die "cannot open input file '$taxanomy' $!";
my $gi=<TAXA>;
chomp $gi;
my @info=split(/\s+/,$gi);
# my $ID;
my $cotu;
my $NOS;
my $Kin;
for (my $i=0; $i<=$#info; $i++)
   {
       if($info[$i] eq "COTU_ID") {$cotu=$i;}
    elsif($info[$i] eq "NO_of_Seqs") {$NOS=$i;}
    elsif($info[$i] eq "Kingdom") {$Kin=$i;}
   }
if(defined $cotu && defined $NOS && defined $Kin) {} else {die "Your input is wrong005!!!\n Please input \"--taxa All_samples_cell_count_annotation_RDP_classifier.txt\"\n $!";}
my %taxa;
while($gi=<TAXA>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        $taxa{$info[$cotu]} = "$info[$Kin];$info[$Kin+1];$info[$Kin+2];$info[$Kin+3];$info[$Kin+4];$info[$Kin+5];NA";
      }
close TAXA;
##Taxanomy finished##

##Bar sequences##
open(FILE, $barseq_data) or die "cannot open input file '$barseq_data' $!";
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
# my $ID;
my $Bar;
my $COTU;
my $seq;
my $link;
my $I1;
my $R2;
for (my $i=0; $i<=$#info; $i++)
   {
       if($info[$i] eq "Bar-sequence-ID") {$Bar=$i;}
    elsif($info[$i] eq "cOTU-ID") {$COTU=$i;}
    elsif($info[$i] eq "Sequence") {$seq=$i;}
    elsif($info[$i] eq "LINK") {$link=$i;}
    elsif($info[$i] eq "I1") {$I1=$i;}
    elsif($info[$i] eq "R2") {$R2=$i;}
   }
if(defined $Bar && defined $COTU && defined $seq && defined $link && defined $I1 && defined $R2) {} else {die "Your input is wrong005!!!\n Please input \"--bar Bar_seq_COTU_clean\"\n $!";}
my $row = 0;
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
while($gi=<FILE>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        my $revcomp = reverse $info[$seq];
        $revcomp =~ tr/ATGCatgc/TACGtacg/;
        print OUTF (">$info[$COTU]($info[$Bar])\t$taxa{$info[$COTU]}\n$revcomp\n");
       $row++;
      }
close FILE;
## Bar sequences##
print "In total there are $row Bar sequences\n";
close OUTF;
print "Done\n";
###Main code###
#####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2023.02.10
