#! /usr/bin/env perl
#################################################################################################
#####Description of this code######
#This code is used to add the cOTU ID to the taxomomies file obtained from BarBIQ_final_tree_rewrite.pl
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_final_Taxonomy_COTU.pl --group file1 --lib file2 --taxa file3
##explaination##
#file1: the output file from BarBIQ_final_Taxonomy_groups.pl
#file2: the output file from BarBIQ_final_lib_COTU_ID.pl
#file3: the output file from BarBIQ_final_tree_rewrite.pl 
########################################################################################################## 
#####Install######
#None
############################################################################################################

#####code#######
use strict;
use warnings;

##save the input file names and output file name##
# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";
##read command##
my ($i,$groupfile, $libary, $taxonomy, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--group")  {$groupfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--lib") {$libary = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--taxa") {$taxonomy = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--group file --lib Libary --taxa taxonomy (--log log_file)\"\n $!";}
    }
if(!$groupfile)   {die "Your input is wrong!!!\n Please input \"--group file\"\n $!";}
if(!$libary)  {die "Your input is wrong!!!\n Please input \"--libary file\"\n $!";}
if(!$taxonomy)  {die "Your input is wrong!!!\n Please input \"--taxa taxonomy\"\n $!";}
##read command##
##check the output name##
my $outputfile=$taxonomy."_COTU";
unlink $outputfile;

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

my %taxa_COTU;
my %no;
open(FILE, $groupfile) or die "cannot open input file '$groupfile' $!";
my $gi=<FILE>;
chomp $gi;
my $g;
my $no;
my $map;
my @info=split(/\s+/,$gi);
if($info[0]=~/\A>/)
    {
      $g=$'; # print "$'\n"; 
      $no=0;
    }
else{die "Your input is wrong!!!004\"\n $!";}
while($gi=<FILE>)
      {
        chomp $gi;
        @info=split(/\s+/,$gi);
        if($info[0]=~/\A>/)
           {
             if($map) {$taxa_COTU{$g}=$map; $no{$g}=$no; $g=$'; undef $map; $no=0;}
            else{die "Your input is wrong!!!005\"\n $!";}
           }
      elsif($info[0] eq "#Taxonomy:")
           {
             $map=join("\t", @info[1..14]);
           }
      else{$no++;}
      }
if($map) {$taxa_COTU{$g}=$map;$no{$g}=$no; undef $map;$no=0;}
else{die "Your input is wrong!!!007\"\n $!";}
close FILE;

my %COTU;
     open(FILE, $libary) or die "cannot open input file '$libary' $!";
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
     if(defined $COTU && defined $seqid) {} else{die "Your input is wrong0015!!!\n";}
     while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          if(!(exists $taxa_COTU{$info[$COTU]}))
             {
               $COTU{$info[$seqid]} = $info[$COTU];
             }
         }
      close FILE; 

open(FILE, $taxonomy) or die "cannot open input file '$taxonomy' $!";
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
$gi=<FILE>;
     chomp $gi;
@info=split(/\s+/,$gi);
my $id;
my $species;
print OUTF ("$gi\tNo_of_Seqs\n");
for(my $i=0; $i<=$#info; $i++)
         {
          if($info[$i] eq "ID") {$id = $i;}
          if($info[$i] eq "species") {$species = $i;}
         }
     if($species) {} else{die "Your input is wrong006!!!\n";}
while($gi=<FILE>)
         {
          chomp $gi;
          @info=split(/\s+/,$gi);
          if(exists $COTU{$info[$id]})
             {
               $info[$species] = $COTU{$info[$id]};
               $info[$id] = $COTU{$info[$id]};
               my $out=join("\t", @info);
               print OUTF ("$out\t1\n");
             }
         }
  foreach my $key (sort keys %taxa_COTU)
        {
          print OUTF ("$key\t$taxa_COTU{$key}\t$no{$key}\n");
        }
  close OUTF;


print "Done\n";
##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.001
##2020.1.7
