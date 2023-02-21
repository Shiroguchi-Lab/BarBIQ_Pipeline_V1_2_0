#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used to correct three types of RepSeqs with indel-related errors: type 1, which has one indel with one substitution; type 2, which has one indel with two substitutions; and type 3, which has two indels.
## The input file is the output file form BarBIQ_sub_chimeras.pl
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_similar.pl --in inputfile --out outputfile
###interpretation##
## --in: the inputfile which is generated from BarBIQ_sub_chimeras.pl.
## --out: outputfile, please set a name for your outputfile.
##########################################################################################################################################
#####Install#####
## please install the perl Modules, IPC::System::Simple, Text::Levenshtein::XS, and Text::WagnerFischer before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system);

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";

##read command##
my $threshold_is=0.2;# a threshold to delete the substitutions+indels
my $threshold_s=0.2; # a threshold to delete the substitutions
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
print "You are using the threshold of $threshold_is to delete the substitutions+indels\n";
print "You are using the threshold of $threshold_s to delete the substitutions\n";

if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
my $outputfile_ratio = "$outputfile"."_ratio";
unlink $outputfile_ratio;
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
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

##Get one cluster's data and find the chimeras##
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!";
chomp($gi=<FILE>);
@info=split(/\s+/,$gi);
my ($cluster_name, $cluster_name_next);
if($info[0] =~ m{\Acluster_([0-9]*)})
  {
   $cluster_name = "cluster_"."$1";
  }
my @sequence;
push @sequence, [@info];
while($gi=<FILE>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    if($info[0] =~ m{\Acluster_([0-9]*)})
        {
          $cluster_name_next = "cluster_"."$1";
        }
    if ($cluster_name_next eq $cluster_name) {
       push @sequence, [@info];
      } else {
       &Find_edit($cluster_name, \@sequence);
       undef @sequence;
       $cluster_name = $cluster_name_next;
       push @sequence, [@info];
      }
    }
&Find_edit($cluster_name, \@sequence);
undef @sequence;
close FILE;

#Find the edit#
sub Find_edit
     {
      my ($cluster, $seq)=@_;
      my ($i,$j, %delete, %mother);
      open (RATIO,'>>', $outputfile_ratio) or die "Could not open file '$outputfile_ratio' $!";
      if ($#{$seq}>=1)
         {
          for($i=0; $i<=($#$seq-1); $i++)
             {
              for($j=$i+1; $j<=$#$seq; $j++)
                 {
                  use Text::Levenshtein::XS qw/distance/;
                  my $distance_I1 = distance($$seq[$i][2],$$seq[$j][2]);
                  my $distance_R2 = distance($$seq[$i][3],$$seq[$j][3]);
                  my $distance=$distance_I1 + $distance_R2;
                  if ($distance <= 5)
                      {
                       use Text::WagnerFischer qw(distance);
                       my $wf_I1 = distance([0,2,1],$$seq[$i][2],$$seq[$j][2]);
                       my $wf_R2 = distance([0,2,1],$$seq[$i][3],$$seq[$j][3]);
                       my $wf=$wf_I1+$wf_R2;
                       my $diff=$wf-$distance;
                       ## Substitutions
                       if ($distance < 5 && $diff == 0)
                          {
                            my $ratio;
                            if($$seq[$i][1] >= $$seq[$j][1])
                               {
                                $ratio=$$seq[$j][1]/$$seq[$i][1];
                                print RATIO ("$$seq[$j][0]\t$$seq[$j][1]\t$$seq[$i][1]\t$ratio\t$distance\t$wf\t$$seq[$j][5]\t$$seq[$i][5]\n");
                                if ($ratio < $threshold_s)
                                   {
                                    if (exists $delete{$$seq[$j][5]})
                                        {print "$$seq[$j][0] has double mothers\n";}
                                   else {$delete{$$seq[$j][5]} = $$seq[$i][5];}
                                    if (exists $mother{$$seq[$i][5]})
                                        {#print "$mother has double sons\n";
                                         $mother{$$seq[$i][5]}=$mother{$$seq[$i][5]}+$$seq[$j][1];
                                        }
                                   else {$mother{$$seq[$i][5]} = $$seq[$j][1];}
                                   }
                               }
                           else{
                                $ratio=$$seq[$i][1]/$$seq[$j][1];
                                print RATIO ("$$seq[$i][0]\t$$seq[$i][1]\t$$seq[$j][1]\t$ratio\t$distance\t$wf\t$$seq[$i][5]\t$$seq[$j][5]\n");
                                if ($ratio < $threshold_s)
                                   {
                                    if (exists $delete{$$seq[$i][5]})
                                        {print "$$seq[$i][0] has double mothers\n";}
                                   else {$delete{$$seq[$i][5]} = $$seq[$j][5];}
                                    if (exists $mother{$$seq[$i][5]})
                                        {#print "$mother has double sons\n";
                                         $mother{$$seq[$j][5]}=$mother{$$seq[$j][5]}+$$seq[$i][1];
                                        }
                                   else {$mother{$$seq[$j][5]} = $$seq[$i][1];}
                                   }
                               }
                           }
                       elsif($distance > 2 && ($diff == 2 || $diff == 4))
                           {
                             my $ratio;
                            if($$seq[$i][1] >= $$seq[$j][1])
                               {
                                $ratio=$$seq[$j][1]/$$seq[$i][1];
                                print RATIO ("$$seq[$j][0]\t$$seq[$j][1]\t$$seq[$i][1]\t$ratio\t$distance\t$wf\t$$seq[$j][5]\t$$seq[$i][5]\n");
                                if ($ratio < $threshold_is)
                                   {
                                    if (exists $delete{$$seq[$j][5]})
                                        {print "$$seq[$j][0] has double mothers\n";}
                                   else {$delete{$$seq[$j][5]} = $$seq[$i][5];}
                                    if (exists $mother{$$seq[$i][5]})
                                        {#print "$mother has double sons\n";
                                         $mother{$$seq[$i][5]}=$mother{$$seq[$i][5]}+$$seq[$j][1];
                                        }
                                   else {$mother{$$seq[$i][5]} = $$seq[$j][1];}
                                   }
                               }
                           else{
                                $ratio=$$seq[$i][1]/$$seq[$j][1];
                                print RATIO ("$$seq[$i][0]\t$$seq[$i][1]\t$$seq[$j][1]\t$ratio\t$distance\t$wf\t$$seq[$i][5]\t$$seq[$j][5]\n");
                                if ($ratio < $threshold_is)
                                   {
                                    if (exists $delete{$$seq[$i][5]})
                                        {print "$$seq[$i][0] has double mothers\n";}
                                   else {$delete{$$seq[$i][5]} = $$seq[$j][5];}
                                    if (exists $mother{$$seq[$i][5]})
                                        {#print "$mother has double sons\n";
                                         $mother{$$seq[$j][5]}=$mother{$$seq[$j][5]}+$$seq[$i][1];
                                        }
                                   else {$mother{$$seq[$j][5]} = $$seq[$i][1];}
                                   }
                               }
                           }
                       ##Substitutions & indels
                      }                      
                 }
             }
         }
      open (OUTF,'>>', $outputfile) or die "Could not open file '$outputfile' $!";
      if(%delete)
          {
           for (my $x=0; $x<=$#{$seq}; $x++)
              {
                if (exists $delete{$$seq[$x][5]}) {}
              else {
                    if (exists $mother{$$seq[$x][5]})
                       {
                        $$seq[$x][1] = $$seq[$x][1] + $mother{$$seq[$x][5]};
                       }
                    my $string = join ("\t", @{$$seq[$x]});
                    print OUTF ("$string\n");
                    }
              }
          }
      else{
           for (my $x=0; $x<=$#{$seq}; $x++)
              {
               my $string = join ("\t", @{$$seq[$x]});
               print OUTF ("$string\n");
              }
          }
      close OUTF;
      close RATIO;
     }
#end of Find the insertion and deletions#

print "Done!!!\n";

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2018.08.28
