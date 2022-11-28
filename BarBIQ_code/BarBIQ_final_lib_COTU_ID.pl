#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to name a unique cell type defined by the determined same-bacterial Bar-sequence(s) a cellular-based Operational Taxonomic Unit (cOTU).
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_lib_COTU_ID.pl --group file1 --lib file2
##explaination##
#file1: output file from BarBIQ_final_overlap_groups.pl
#file2: the output file from BarBIQ_final_libary.pl
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
my ($i,$groupfile, $libary, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--group")  {$groupfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--bar") {$libary = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];} 
     else                         {die "Your input is wrong!!!\n Please input \"--group file --bar bar-sequence-file (--log log_file)\"\n $!";}
    }
if(!$groupfile)   {die "Your input is wrong!!!\n Please input \"--group file\"\n $!";}
if(!$libary)  {die "Your input is wrong!!!\n Please input \"--bar bar-sequence-file\"\n $!";}
##read command##
##check the output name##
my $outputfile=$libary."_COTU";
unlink $outputfile;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

open(FILE, $libary) or die "cannot open input file '$libary' $!";
  my $gi=<FILE>;
  $gi=<FILE>;
  chomp $gi;
  my @info=split(/\s+/,$gi);
  my $dataset;
  if($info[0] =~ /\ABar-sequence-/) {$dataset= "";} else{ die "Your input is wrong008!!!\n";}
  close FILE;

my %group;
     open(FILE, $groupfile) or die "cannot open input file '$groupfile' $!";
     while(my $gi=<FILE>)
         {
           chomp $gi;
           my @info=split(/\s+/,$gi);
           for(my $i=1; $i<=$#info; $i++)
              {
                $group{$info[$i]} = $info[0];
              }
         }
     close FILE;

     open(FILE, $libary) or die "cannot open input file '$libary' $!";
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
          if(exists $group{$info[$seqid]}) 
                 { 
                   if($group{$info[$seqid]} =~ /\AGroup_/)
                       {  
                          my $ID=2000+$';
                          $info[$COTU] = "COTU-$ID"; 
                          my $out=join("\t", @info); print OUTF ("$out\n");
                        }
                     else{ die "Your input is wrong006!!!\n";}
                 }
         else{ 
               if($info[$seqid] =~ /\ABar-sequence/)  
                       {  $info[$COTU] = "COTU$'"; my $out=join("\t", @info); print OUTF ("$out\n");}
                    else{ die "Your input is wrong007!!!\n";}
             }
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
