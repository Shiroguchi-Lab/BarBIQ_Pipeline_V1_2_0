#! /usr/bin/env perl
###############################################################
#######Description of this code#####
#This code is used to write the information including raw counts, cellular concentration, taxonomies, etc. into a excel file.
##########################################################################
######how to run this code #####
##command##
#BarBIQ_final_bac_con_nomalization_xlsx.pl --count file1 --total file2 --out file3 --taxa method --dataname name
##explaination##
#file1: the output file form BarBIQ_M_final_add_RDP_prediction_Taxonomy_to_bac_count.pl or BarBIQ_final_add_mapping_results_bac_counting.pl, which the adding taxonomy average bacterial counting file
#file2: Total concentration of each sample, should be prepared according the example (example_ Total_concentration.txt) we provided.
#file3: outputfile.xlsx
#method: RDP/GG/Silva/Classifier 
#name: a dataname used to label the cOTU ID, e.g., “MK” in COTU-MK-01
########################################################################################################
#####Install#####
#Please install the module Excel::Writer::XLSX before using this code
##############################################################################################

#####code#####

use strict;
use warnings;
use Excel::Writer::XLSX;

##read command##
my ($i,$ave_data,$totalcon,$raw_data, $index, $outputfile,$taxa, $dtnm, $cutoff, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--count") {$ave_data = $ARGV[$i+1];}
   #  elsif ($ARGV[$i] eq "--raw") {$raw_data = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--total") {$totalcon = $ARGV[$i+1];}
#     elsif ($ARGV[$i] eq "--index") {$index = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--taxa") {$taxa = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--cutoff") {$cutoff = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
#     elsif ($ARGV[$i] eq "--dataname") {$dtnm = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--count cell_ccount_data --total total_concentration_file --out outputfile --taxa RDP/GG/Silva/Classifier (--log log_file)\"\n $!";}   
    }
if(!$totalcon)   {die "Your input is wrong!!!\n Please input \"--total total_concentration_file\"\n $!";}
if(!$ave_data)  {die "Your input is wrong!!!\n Please input \"--count cell_ccount_data\"\n $!";}
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
###Main code###
open(TOTAL, $totalcon) or die "cannot open input file '$totalcon' $!";
my $gi=<TOTAL>;
chomp $gi;
my $totcon;
my $sample_name;
my @info=split(/\s+/,$gi);
for (my $i=0; $i<=$#info; $i++)
    {
      if ($info[$i] eq "Sample") {$sample_name=$i;}
      if ($info[$i] eq "Conc.(copies/mg)") {$totcon = $i;}
    }
if ((! defined $sample_name) || (! defined $totcon)) {die "Your input $totalcon is wrong!!!002\"\n $!";}
my %totolconc;
while($gi=<TOTAL>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        $totolconc{$info[$sample_name]} = $info[$totcon];
      }
close TOTAL;

## print the introduction
my $workbook  = Excel::Writer::XLSX->new( $outputfile );
my $worksheet = $workbook->add_worksheet('ReadMe');
$worksheet->write(0, 0, "Title for sheets");
$worksheet->write(1, 0, "Sequencing-determined counts:");
$worksheet->write(1, 1, "Cell numbers of each cOTU determined by sequenceing");
$worksheet->write(2, 0, "Proportional-concentration:");
$worksheet->write(2, 1, "Proportional concentration of each cOTU obtained from the cell numbers of each cOTU determined by sequencing divided by the total number of cells in the sample");
$worksheet->write(3, 0, "Total-concentration:");
$worksheet->write(3, 1, "Total concentration of each sample measured by ddPCR");
$worksheet->write(4, 0, "Absolute-concentration:");
$worksheet->write(4, 1, "Absolute cell concentration of each cOTU obtained from the cell numbers of each cOTU determined by sequencing normalized by the total concentration of the sample.");
$worksheet->write(6, 0, "Legend");
$worksheet->write(7, 0, "Kingdom/Phylum/Class/Order/Family/Genus:");
if($taxa eq "Silva")
{
$worksheet->write(7, 1, "The Taxonomy: mapped to the Silva database; \(* - - -\) means it cannot be determined from this level, because it can be mapped to different names from this level; If \(*\) at Kingdom level, it means it cannot be mapped to any references");
}
elsif($taxa eq "RDP")
{
$worksheet->write(7, 1, "The Taxonomy: mapped to the RDP database; \(* - - -\) means it cannot be determined from this level, because it can be mapped to different names from this level; If \(*\) at Kingdom level, it means it cannot be mapped to any references");
}
elsif($taxa eq "GG")
{
$worksheet->write(7, 1, "The Taxonomy: mapped to the GREENGENE database; \(* - - -\) means it cannot be determined from this level, because it can be mapped to different names from this level; If \(*\) at Kingdom level, it means it cannot be mapped to any references");
}
elsif($taxa eq "Classifier")
{
$worksheet->write(7, 1, "Taxonomies of each cOTU; the taxonomies of each Bar sequence in the given cOTU were predicted by RDP Classifier; the highest scored prediction (score shown in parentheses) of all Bar sequences in the given cOTU was chosen to be the prediction of the cOTU");
}
else{die "Your input is wrong!!!\n Please input \"--taxa RDP/GG/Silva/Classifier\"\n $!";}
$worksheet->write(8, 0, "cOTU-ID:");
$worksheet->write(8, 1, "ID of cOTU");
$worksheet->write(9, 0, "No_of_Bar-sequence:");
$worksheet->write(9, 1, "Number of Bar sequences in the given cOTU");
# $worksheet->write(7, 0, "Raw_data: Raw counts from BarBIQ of all technique repeats of all samples;");
# $worksheet->write(8, 0, "Index_IF: the information about the index of raw data to the samples;");

open(FILE, $ave_data) or die "cannot open input file '$ave_data' $!";
$worksheet = $workbook->add_worksheet('Sequencing-determined counts');

$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
# my $ID;
my $Kingdom;
my $No_seqs;
my $COTU;
my %total_number;
my @samples;
for (my $i=0; $i<=$#info; $i++)
   {
    if($info[$i] eq "NO_of_Seqs") {$No_seqs=$i;}
    if($info[$i] eq "Kingdom") {$Kingdom=$i;}
    if($info[$i] eq "COTU_ID") {$COTU=$i;}
   }
if (! defined $Kingdom ) {die "Your input is wrong!!!003\"\n $!";}
if (! defined $COTU) {die "Your input is wrong!!!004\"\n $!";}
if (! defined $No_seqs) {die "Your input is wrong!!!005\"\n $!";}
print "You have samples:\n";
for (my $i=1; $i<$No_seqs; $i++)
   {
     $total_number{$info[$i]} = 0;
     $samples[$i] = $info[$i];
     print "$info[$i]\n";
   }
$worksheet->write(0, 0, "cOTU_ID");
for (my $i=1; $i<$No_seqs; $i++)
   {
      $worksheet->write(0, $i, $info[$i]);
   }
my $row = 1;
while($gi=<FILE>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        if ($info[0] =~ /\ACOTU/ )
           {
             my $ids = $info[0]; # "COTU-".$dtnm."-".$';
             $worksheet->write($row, 0, $ids);
           }
        else
           {die "Your input is wrong!!!006\"\n $!";}
        for (my $i=1; $i<$No_seqs; $i++)
          {   
            $worksheet->write($row, $i, $info[$i]);
          }
       for (my $i=1; $i<$No_seqs; $i++)
          {
            $total_number{$samples[$i]} = $total_number{$samples[$i]} + $info[$i];
          }
       $row++;
      }
close FILE;
## Proportinal concentration

$worksheet = $workbook->add_worksheet('Proportinal-concentration');
open(FILE, $ave_data) or die "cannot open input file $ave_data' $!";
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
$worksheet->write(0, 0, "cOTU_ID");
for (my $i=1; $i<$No_seqs; $i++)
   {
     $worksheet->write(0, $i, $info[$i]);
   }
$row = 1;
while($gi=<FILE>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        if ($info[0] =~ /\ACOTU/ )
           {
             my $ids = $info[0]; # "COTU-".$dtnm."-".$';
             $worksheet->write($row, 0, $ids);
           }
        else        
           {die "Your input is wrong!!!006\"\n $!";}
        for (my $i=1; $i<$No_seqs; $i++)
          {
            my $conc=$info[$i]/$total_number{$samples[$i]};
            $worksheet->write($row, $i, $conc);
          }
      # for (my $i=$No_seqs; $i<=$#info; $i++)   {$worksheet->write($row, $i, $info[$i]); }
       $row++;
      }
close FILE;

##Total conc
$worksheet = $workbook->add_worksheet('Total-concentration');
open(TOTAL, $totalcon) or die "cannot open input file '$totalcon' $!";
$row = 0;
while($gi=<TOTAL>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        for (my $i=0; $i<=$#info; $i++)
           {
             $worksheet->write($row, $i, $info[$i]);
           }
        
        $row++;
      }
close TOTAL;

## Absolute concentration
my @taxa = ("Kingdom","Phylum",	"Class", "Order", "Family", "Genus");
$worksheet = $workbook->add_worksheet('Absolute-concentration');
open(FILE, $ave_data) or die "cannot open input file '$ave_data' $!";
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
$worksheet->write(0, 0, "cOTU_ID");
for (my $i=1; $i<$No_seqs; $i++)
   {
     $worksheet->write(0, $i, $info[$i]);
   }
$worksheet->write(0, $No_seqs, "No_of_BISes");
for (my $i=$No_seqs+1; $i<=$#info; $i++)
   {
     $worksheet->write(0, $i, $info[$i]);
   }
$worksheet->write(0, $No_seqs+7, "Lowest_level_taxon_with_score>=$cutoff");
$row = 1;
while($gi=<FILE>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
         if ($info[0] =~ /\ACOTU/ )
           {
             my $ids = $info[0]; # "COTU-".$dtnm."-".$';
             $worksheet->write($row, 0, $ids);
           }
        else
           {die "Your input is wrong!!!006\"\n $!";}
        for (my $i=1; $i<$No_seqs; $i++)
          {
          #   print "$info[$i]\n$total_number{$samples[$i]}\n$totolconc{$samples[$i]}\n";
            my $conc=$info[$i]/$total_number{$samples[$i]}*$totolconc{$samples[$i]};
            $worksheet->write($row, $i, $conc);
      #      last;
          }
       $worksheet->write(0, $No_seqs, $info[$No_seqs]);
       my $low_taxon = "Bacteria";
       my $low_level = "Kingdom";
       for (my $i=$No_seqs+1; $i<=$#info; $i++)
          {
            if($info[$i] =~ /(.*)\((.*)\)/)
                {
                  if($2 >= $cutoff) { $low_taxon = $info[$i]; $low_level = $taxa[$i-$No_seqs-1];}
                }
             else{die "Your input is wrong!!!009\"\n $!";}
            $worksheet->write($row, $i, $info[$i]);
          }
       $worksheet->write($row, $No_seqs+7, "$low_level"."_$low_taxon");
       $row++;
      }
close FILE;

##Raw data
# $worksheet = $workbook->add_worksheet('Raw_data');
# open(RAW, $raw_data) or die "cannot open input file '$raw_data' $!";
# $gi=<RAW>;
# chomp $gi;
# @info=split(/\s+/,$gi);
# my $ID;
# if ($info[0] eq "ID" && $info[$#info] =~ /\AS/) {} else {die "your input file $raw_data is wrong0019\n";}
# print "You have indexes:\n";
# for (my $i=0; $i<=$#info; $i++)
#  {
#   $worksheet->write(0, $i, $info[$i]);
#   print "$info[$i]\n";
#  }
#$row = 1;
#while($gi=<RAW>)
#      {
#        chomp $gi;
#        @info=split(/\s+/,$gi);
#        for (my $i=0; $i<=$#info; $i++)
#          {
#            $worksheet->write($row, $i, $info[$i]);
#          }
#       $row++;
#      }
# close RAW;

##Index file
# $worksheet = $workbook->add_worksheet('Index_IF');
# open(INDEX, $index) or die "cannot open input file '$index' $!";
# $row = 0;
# while($gi=<INDEX>)
#      {
#        chomp $gi;
#        @info=split(/\s+/,$gi);
#        for (my $i=0; $i<=$#info; $i++)
#           {
#             $worksheet->write($row, $i, $info[$i]);
#           }
#        $row++;
#      }
# close INDEX;

$workbook->close;
print "Done\n";
###Main code###
#####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2019.01.05
