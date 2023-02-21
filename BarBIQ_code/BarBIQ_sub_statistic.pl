#! /usr/bin/env perl
######Description of this code#####
## This code is used for each unique RepSeq (RepSeq type), to count the number of BClusters containing the given RepSeq type in each index.
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_statistic.pl --in inputfile --out outputfile
###explaination##
## --in: the inputfile is the output file from BarBIQ_sub_similar.pl.
## --out: outputfile, please set a name for your outputfile
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
my ($i,$inputfile,$outputfile, $log_file);
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --out outputfile (--log log_file_name)\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out: outputfile\"\n $!";}

if(!$log_file) {
                 $log_file = $outputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
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
print "Inputfiles are OK!\nStart to calculating:\n";
##check the inputfile##

###Main code###
print "start to statistic the inputfile......\n";
my %inputfile_stat;
my %inputfile_Link;
open (FILER,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
while($gi=<FILER>)
   {
     chomp $gi;
     @info=split(/\s+/, $gi);
     if(exists $inputfile_stat{$info[5]})
        {
          $inputfile_stat{$info[5]}++;
        }
    else{
          $inputfile_stat{$info[5]}=1;
          my $length_I1=length($info[2]);
          my $length_R2=length($info[3]);
          $inputfile_Link{$info[5]}="$info[4]\t$length_I1\t$length_R2";
        }
   }
close FILER;
open(OUTFS, '>>', $outputfile) or die "canot open input file '$outputfile' $!";
my @stat;
foreach my $key(keys %inputfile_stat)
   {
     push @stat,[$inputfile_stat{$key},$key,$inputfile_Link{$key}];
   }

my @out = sort {$a->[0] <=> $b->[0]} @stat;
$i=1;
for my $out(@out)
  {
  my $string=join("\t",@$out);
  print OUTFS ("$i\t$string\n");
  $i++;
  }
close OUTFS;
#Review the original sequence and keep all right sequences##

print "Done!!!\n";
###Main code###
####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
