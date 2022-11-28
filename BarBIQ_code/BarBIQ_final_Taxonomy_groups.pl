#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to predict the taxonomies for the cOTUs which contain mutiple Bar sequences based on the predicted taxonomies of each Bar sequence using RDP classifier.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_Taxonomy_groups.pl --group file1 --taxa file2
##explaination##
#file1: the output file from BarBIQ_final_overlap_groups.pl
#file2: the output file from BarBIQ_final_tree_rewrite.pl
####################################################################################################

#####code#######
use strict;
use warnings;

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##read command##
my ($i,$groupfile, $taxonomy, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--group")  {$groupfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--taxa") {$taxonomy = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}     
     else                         {die "Your input is wrong!!!\n Please input \"--group file --taxa taxonomy (--log log_file)\"\n $!";}
    }
if(!$groupfile)   {die "Your input is wrong!!!\n Please input \"--group file\"\n $!";}
if(!$taxonomy)  {die "Your input is wrong!!!\n Please input \"--taxa taxonomy\"\n $!";}
##read command##
##check the output name##
my $outputfile="$groupfile"."_RDPpdt_taxonomy";
unlink $outputfile;
if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";


my %mapping_bac;
     open(FILE, $taxonomy) or die "cannot open input file '$taxonomy' $!";
     my $gi=<FILE>;
     chomp $gi;
     my @info=split(/\s+/,$gi);
     my $bac;
     my $seqid;
     for(my $i=0; $i<=$#info; $i++)
         {
          if($info[$i] eq "domain") {$bac = $i;}
          if($info[$i] eq "ID") {$seqid = $i;}
         }
     if($bac) {} else{die "Your input is wrong001!!!\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          if(exists $mapping_bac{$info[$seqid]}) {print "$gi\n";}
         else{ $mapping_bac{$info[$seqid]} = "$info[$seqid]\t$info[$bac]\t$info[$bac+1]\t$info[$bac+2]\t$info[$bac+3]\t$info[$bac+4]\t$info[$bac+5]\t$info[$bac+6]\t$info[$bac+7]\t$info[$bac+8]\t$info[$bac+9]\t$info[$bac+10]\t$info[$bac+11]";}
         }
      close FILE; 
     
use Data::Dumper;
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
open(FILE1, $groupfile) or die "cannot open input file '$groupfile' $!";
while($gi=<FILE1>)
           {
             chomp $gi;
             @info=split(/\s+/,$gi);
             print OUTF (">$info[0]\n");
             my @common;
             for(my $i=1; $i<=$#info; $i++)
                {
                  if (exists $mapping_bac{$info[$i]})
                     {
                       print OUTF ("$mapping_bac{$info[$i]}\n");
                       my @xxx=split(/\s+/,$mapping_bac{$info[$i]});
                       if(@common)
                           {
                             for(my $c=0; $c<=$#common; $c=$c+2)
                                {
                                  if($common[$c] ne $xxx[$c+1])
                                     {
                                        if($common[$c+1] < $xxx[$c+2])
                                             { 
                                                @common = @xxx[1..12]; last;
                                             }
                                     }
                                 else{
                                       if($common[$c+1] < $xxx[$c+2] && $c == 10) { @common = @xxx[1..12];}
                                     }
                                }
                           }
                       else{
                             @common = @xxx[1..12];
                           }
                     } 
                else { print ("not exist in the lib:\n $gi\n");}
                }
             if(@common) {
             my $common=join("\t", @common);
             print OUTF ("#Taxonomy:\t$common\t$info[0]\t0\n");
            } else { print OUTF ("#Taxonomy:\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n");}
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
