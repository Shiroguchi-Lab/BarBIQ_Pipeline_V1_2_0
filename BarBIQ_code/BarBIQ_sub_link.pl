#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is use to link I1 and R2 RepSeqs by their overlapping 3’ end sequences, if they are linked, labeled as "LINK", if not, labeled as "I1R2".
## they will be linked when the number of overlapped bases was more than 5 bases;
## the inputfile is the output file from BarBIQ_sub_shift.pl
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_link.pl --in inputfile --out outputfile
###explaination##
## --in: the inputfile which is generated from BarBIQ_sub_shift.pl.
## --out: outputfile, please set a name for your outputfile
##########################################################################################################################################
#####Install#####
## None
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
my $outputfile_2="$outputfile"."_bases"; ### save how many bases are overlapped for each rep-seq
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfile is:\n$inputfile\n";
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
if(!(($#info == 3) && ($info[0] =~ /\Acluster_/)))
   {
    die "Your input file $inputfile is wrong!!! please check!!!\n $!";
   }
if ($info[3] eq "No")
   {
    print "This is a single end data!!!please use another code!!!\n";
    die;
   }
close FILE;
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

##Find the overlap##
my $overlapcon=5; # how many bases should be considered as overlapped.
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!";
open (OUTF,'>>', $outputfile) or die "Could not open file '$outputfile' $!";
open (OUTF2,'>>', $outputfile_2) or die "Could not open file '$outputfile_2' $!";
while($gi=<FILE>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    my ($overlap, $seq, $bases)=&Find_overlap($info[2], $info[3]);
    print OUTF ("$gi\t$overlap\t$seq\n");    
    print OUTF2 ("$bases\t$info[2]\t$info[3]\n");
    }
close FILE;
close OUTF;
close OUTF2;

#find the overlapped sequence#
  sub Find_overlap 
      {
       my ($I1_seq, $R2_seq)=@_;
       my @seq_I1 = split("", $I1_seq);
       my @seq_R2 = split("", $R2_seq);
       my @seq_R2_rev;
       my $position=0;
       for($i=0; $i<=$#seq_R2; $i++)
        {
            if ($seq_R2[$i] eq "A") {$seq_R2_rev[$#seq_R2-$i]="T";}
         elsif ($seq_R2[$i] eq "T") {$seq_R2_rev[$#seq_R2-$i]="A";}
         elsif ($seq_R2[$i] eq "C") {$seq_R2_rev[$#seq_R2-$i]="G";}
         elsif ($seq_R2[$i] eq "G") {$seq_R2_rev[$#seq_R2-$i]="C";}
         else {print "Something is wrong!!!\n"}
        }
       my $I1_length=$#seq_I1;
       my $R2_length=$#seq_R2;
       my $start;
       if ($I1_length>=$R2_length) {$start=$I1_length-$R2_length;}
       else {$start=0;}
       for($i=$start; $i<$#seq_I1-$overlapcon; $i++)
           {
            if($seq_I1[$i] eq $seq_R2_rev[0])
                {
                 my @seq1= @seq_I1;
                 my @end_seq = splice @seq1, $i;
                 my $end_seq = join("", @end_seq);
                 my @seq2=@seq_R2_rev;
                 splice @seq2, $#end_seq+1;
                 my $head_seq = join("", @seq2);
                 if ($end_seq eq $head_seq)
                     {
                      $position = $i;
                      last;
                     }
                } 
             if($position !=0)
                 {last;}
           }
       my ($seq_I1,$seq_R2,$seq_I1R2);
       if ($position != 0)
            {
             splice @seq_I1, $position;
             $seq_I1 = join("", @seq_I1);
             $seq_R2 = join("", @seq_R2_rev);
             $seq_I1R2="$seq_I1"."$seq_R2";
             my $base=length($seq_I1)+length($seq_R2)-length($seq_I1R2);
             return ("LINK", $seq_I1R2, "$base");
            }
       else
           {
            $seq_I1 = join("", @seq_I1);
            $seq_R2 = join("", @seq_R2_rev);
            $seq_I1R2="$seq_I1"."$seq_R2";
            return ("I1R2", $seq_I1R2, "0");
           }
      }
#find the overlapped sequence#
print "Done!!!\n";

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
