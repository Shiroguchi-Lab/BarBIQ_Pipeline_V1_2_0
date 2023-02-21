#! /usr/bin/env perl
###############################################################
#######Description of this code#####
#This code is used to add predicted taxonomies by RDP classifier to the final cOTU counts file.
##########################################################################
######how to run this code ######
##command##
#BarBIQ_final_add_Taxonomy_to_bac_count_publish.pl --bacc file1 --taxa file2
###explaination##
#file1: the out put file from BarBIQ_final_compare_datasets.pl 
#file2: the out put file from BarBIQ_final_Taxonomy_COTU.pl
#########################################################################################
#####Install######
#None
##########################################################################################

#####code#####

use strict;
use warnings;

##read command##
my ($i,$bac_count,$taxanomy,$groupfile, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--bacc")  {$bac_count = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--taxa")  {$taxanomy = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log")  {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--bacc bac_count --taxa taxanomy (--log log_file)\"\n $!";}
    }
if(!$bac_count)   {die "Your input is wrong!!!\n Please input \"--bacc bac_count\"\n $!";}
if(!$taxanomy)   {die "Your input is wrong!!!\n Please input \"--taxa taxanomy\"\n $!";}
my $outputfile = $bac_count."_annotation_RDP_classifier.txt";
unlink $outputfile;
if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
#
#
###Main code###
my %taxa_bac;
my %No_seqs;
     open(FILE, $taxanomy) or die "cannot open input file '$taxanomy' $!";
     my $gi=<FILE>;
     chomp $gi;
     my $title=$gi;
     my @info=split(/\s+/,$gi);
     my $bac;
     my $seqid;
     my $No_seq;
     for(my $i=0; $i<=$#info; $i++)
         {
          if($info[$i] eq "domain") {$bac = $i;}
          if($info[$i] eq "ID") {$seqid = $i;}
          if($info[$i] eq "No_of_Seqs") {$No_seq = $i;}
         }
     if($bac) {} else{die "Your input is wrong001!!!\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          if(exists $taxa_bac{$info[$seqid]}) {print "$gi\n";}
         else{ 
               my @taxas;
               for(my $i=1; $i<=$No_seq-3; $i=$i+2)
                  {
                    push @taxas, "$info[$i]($info[$i+1])";
                  }
               $taxa_bac{$info[$seqid]} = join("\t", @taxas);
               $No_seqs{$info[$seqid]} = $info[$No_seq];
             }
         }
      close FILE;

open(FILE, $bac_count) or die "cannot open input file '$bac_count' $!";
open(OUTF, '>>', $outputfile)  or die "canot open input file '$outputfile' $!";
$gi=<FILE>; chomp $gi;
@info=split(/\s+/,$gi);
shift @info;
my $title2=join("\t", @info);
print OUTF ("COTU_ID\t$title2\tNO_of_Seqs\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\n");
while($gi=<FILE>)
{
     chomp $gi;
     @info=split(/\s+/,$gi);
     my $ID=shift @info;
     if (exists $taxa_bac{$ID})
        {
#         print "$info[5]\n";
         my $data=join("\t", @info);
         print OUTF ("$ID\t$data\t$No_seqs{$ID}\t$taxa_bac{$ID}\n");
        }
    else{
         print ("$ID is not existed in the taxonomy file\n");
        }
}
close FILE;
close OUTF;
print "Done\n";
###Main code###
#####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2018.11.20

