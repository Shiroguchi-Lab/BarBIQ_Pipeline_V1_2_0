#! /usr/bin/env perl
###############################################################
#######Description of this code#####
#This code is used to write the Bar sequences into a excel file.
##########################################################################
######how to run this code #####
##command##
#BarBIQ_final_bar_seq_xlsx.pl --bar file1 --out file1
##explaination##
#file1: the output file from  
#file2: outputfile.xlsx
########################################################################################################
#####Install#####
#Please install the module Excel::Writer::XLSX before using this code
##############################################################################################

#####code#####

use strict;
use warnings;
use Excel::Writer::XLSX;

##read command##
my ($i,$barseq_data,$outputfile,$log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--bar") {$barseq_data = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
#     elsif ($ARGV[$i] eq "--dataname") {$dtnm = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--bar Bar_seq_COTU_clean --out outputfile (--log log_file)\"\n $!";}   
    }
if(!$barseq_data)  {die "Your input is wrong!!!\n Please input \"--bar Bar_seq_COTU_clean\"\n $!";}
# if(!$dtnm)  {die "Your input is wrong!!!\n Please input \"--dataname XX\"\n $!";}
# if(!$raw_data)  {die "Your input is wrong!!!\n Please input \"--raw raw_file\"\n $!";}
# if(!$index)  {die "Your input is wrong!!!\n Please input \"--index index_file\"\n $!";}
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

##Read me##
my $workbook  = Excel::Writer::XLSX->new( $outputfile );
my $worksheet = $workbook->add_worksheet('ReadMe');
$worksheet->write(0, 0, "BarBIQ-identified sequences (Bar sequences)");
$worksheet->write(1, 0, "Legend");
$worksheet->write(2, 0, "Bar-sequence-ID:");
$worksheet->write(2, 1, "ID of Bar sequence");
$worksheet->write(3, 0, "cOTU-ID:");
$worksheet->write(3, 1, "ID of cOTU");
$worksheet->write(4, 0, "Sequence(5\'->3\'):");
$worksheet->write(4, 1, "Seqeunce of Bar sequence containning both I1 and R2 and linked by their overlopped sequences");
$worksheet->write(5, 0, "LINK");
$worksheet->write(5, 1, "\"LINK\" indicates I1 and R2 are overlapped; \"I1R2\" indicates I1 and R2 are not overlapped");
$worksheet->write(6, 0, "Sequence-I1(5\'->3\')");
$worksheet->write(6, 1, "Seqeunce of Bar sequence identified by I1 reads");
$worksheet->write(7, 0, "Sequence-R2(5\'->3\')");
$worksheet->write(7, 1, "Seqeunce of Bar sequence identified by R2 reads");
##Read me##

##Bar sequences##
open(FILE, $barseq_data) or die "cannot open input file '$barseq_data' $!";
$worksheet = $workbook->add_worksheet('Bar-sequences');

my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
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
$worksheet->write(0, 0, "Bar-sequence-ID");
$worksheet->write(0, 1, "cOTU-ID");
$worksheet->write(0, 2, "Sequence(5\'->3\')");
$worksheet->write(0, 3, "LINK");
$worksheet->write(0, 4, "Sequence-I1(5\'->3\')");
$worksheet->write(0, 5, "Sequence-R2(5\'->3\'");
my $row = 0;
while($gi=<FILE>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        $worksheet->write($row, 0, $info[$Bar]);
        $worksheet->write($row, 1, $info[$COTU]);
        my $revcomp = reverse $info[$seq];
        $revcomp =~ tr/ATGCatgc/TACGtacg/;
        $worksheet->write($row, 2, $revcomp);
        $worksheet->write($row, 3, $info[$link]);
        my $I1_seq = substr($info[$seq], 0, $info[$I1]);
        $I1_seq = reverse $I1_seq;
        $I1_seq =~ tr/ATGCatgc/TACGtacg/;
        $worksheet->write($row, 4, $I1_seq );
        my $R2_seq = substr($info[$seq], (length($info[$seq])-$info[$R2]), $info[$R2]);
        $R2_seq = reverse $R2_seq;
        $R2_seq =~ tr/ATGCatgc/TACGtacg/;
        $worksheet->write($row, 5, $R2_seq);  
       $row++;
      }
close FILE;
## Bar sequences##
print "In total there are $row Bar sequences\n";

$workbook->close;
print "Done\n";
###Main code###
#####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2022.11.25
