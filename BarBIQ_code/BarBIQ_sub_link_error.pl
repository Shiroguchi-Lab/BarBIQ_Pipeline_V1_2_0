#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used to remove one substitution errors.
## The input file is the output file from BarBIQ_sub_link.pl
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_link_error.pl --in inputfile --out outputfile
###interpretation##
## --in: the inputfile which is generated from BarBIQ_sub_link.pl.
## --out: outputfile, please set a name for your outputfile
##########################################################################################################################################
#####Install#####
## please install the perl Module Text::Levenshtein::XS before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use Text::Levenshtein::XS qw/distance/;

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
my $outputfile_2="$outputfile"."_errors"; ### save how many bases are overlapped for each rep-seq
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfile is:\n$inputfile\n";
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
if(!(($#info == 5) && ($info[0] =~ /\Acluster_/)))
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

##Find the errors##
my %possible_error;
my %if;
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
       &Find_error($cluster_name, \@sequence);
       undef @sequence;
       $cluster_name = $cluster_name_next;
       push @sequence, [@info];
      }
    }
&Find_error($cluster_name, \@sequence);
undef @sequence;
close FILE;

### review all sequences
my (%pair_b, %pair_s, %pair_n);
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!";
chomp($gi=<FILE>);
@info=split(/\s+/,$gi);
if($info[0] =~ m{\Acluster_([0-9]*)})
  {
   $cluster_name = "cluster_"."$1";
  }
my %sequence;
$sequence{$info[5]} = $info[0];
while($gi=<FILE>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    if($info[0] =~ m{\Acluster_([0-9]*)})
        {
          $cluster_name_next = "cluster_"."$1";
        }
    if ($cluster_name_next eq $cluster_name) {
       $sequence{$info[5]} = $info[1];
      } else {
       &Review_error(\%sequence);
       undef %sequence;
       $cluster_name = $cluster_name_next;
       $sequence{$info[5]} = $info[1];
      }
    }
&Review_error(\%sequence);
undef %sequence;
close FILE;

#### delete the junk
my %delete;
open (OUTF2,'>>', $outputfile_2) or die "Could not open file '$outputfile_2' $!";
foreach my $key (keys %possible_error)
    {
       print OUTF2 ("$pair_b{$key}\t$pair_s{$key}\t$pair_n{$key}\t$key\t$possible_error{$key}\n");
       if($pair_b{$key} > 0)
          {
            $delete{$key}= $possible_error{$key};
          }
    }
my @no=keys %delete;
# print "$#no\n";
close OUTF2;

open (FILE,$inputfile) or die "Could not open file '$inputfile' $!";
open (OUTF,'>>', $outputfile) or die "Could not open file '$outputfile' $!";
chomp($gi=<FILE>);
@info=split(/\s+/,$gi);
if($info[0] =~ m{\Acluster_([0-9]*)})
  {
   $cluster_name = "cluster_"."$1";
  }
undef %sequence;
$sequence{$info[5]} = $gi;
while($gi=<FILE>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    if($info[0] =~ m{\Acluster_([0-9]*)})
        {
          $cluster_name_next = "cluster_"."$1";
        }
    if ($cluster_name_next eq $cluster_name) {
       $sequence{$info[5]} = $gi;
      } else {
       &Delete($cluster_name, \%sequence);
       undef %sequence;
       $cluster_name = $cluster_name_next;
       $sequence{$info[5]} = $gi;
      }
    }
&Delete($cluster_name, \%sequence);
undef %sequence;
close FILE;
close OUTF;

##delete
sub Delete
    {
      my ($csname, $seqseq)=@_;
      foreach my $key (keys %{$seqseq})
        {
          if(exists $delete{$key})
             {
              if(exists $$seqseq{$delete{$key}})
                 {
                   my @d=split(/\s+/,$$seqseq{$key});
                   my @m=split(/\s+/,$$seqseq{$delete{$key}});
                   $m[1] = $m[1] + $d[1];
                   $$seqseq{$delete{$key}} = join("\t", @m);
                   delete $$seqseq{$key};
                  }
              else{
                    my @d=split(/\s+/,$$seqseq{$key});
                    my $dd="$d[0]\t$d[1]\t$if{$delete{$key}}";
                    delete $$seqseq{$key};
                    $$seqseq{$delete{$key}} = $dd;
                  }
             }
        }
     foreach my $key (keys %{$seqseq})
       {
         print OUTF ("$$seqseq{$key}\n");
       }
    }
#review
 sub Review_error
    {
      my ($seqseqseq)=@_;
      my %consider;
      foreach my $key (keys %{$seqseqseq})
         {
           if(exists $possible_error{$key}) { $consider{$key} = $possible_error{$key}; }
         }
     if(%consider)
        {
          foreach my $key (keys %consider)
           {
            if(exists $$seqseqseq{$consider{$key}}) 
                {
                  if( $$seqseqseq{$consider{$key}} >= $$seqseqseq{$key} ) 
                       {
                        if(exists $pair_b{$key})
                            {
                             $pair_b{$key}++;
                            }
                       else{
                             $pair_b{$key}=1;
                             $pair_s{$key}=0;
                             $pair_n{$key}=0;
                           }
                       }
                   else{
                          if(exists $pair_b{$key})
                            {
                             $pair_s{$key}++;
                            }
                       else{
                             $pair_b{$key}=0;
                             $pair_s{$key}=1;
                             $pair_n{$key}=0;
                           }
                       }
                    }
                else{
                      if(exists $pair_b{$key})
                            {
                             $pair_n{$key}++;
                            }
                       else{
                             $pair_b{$key}=0;
                             $pair_s{$key}=0;
                             $pair_n{$key}=1;
                           }
                    }
          }
        }
    }
#find the error overlapped sequence#
  sub Find_error 
      {
       my ($csname, $seqseq)=@_;
       my $need1=0;
       my $need2=0;
       for(my $i=0; $i<=$#$seqseq; $i++)
          {
            if($$seqseq[$i][4] eq "I1R2") { $need1=1;}
            if($$seqseq[$i][4] eq "LINK") { $need2=1;}
          }
       if($need1 && $need2)
          {
            for(my $i=0; $i<=$#$seqseq; $i++)
               {
                 for(my $j=0; $j<=$#$seqseq; $j++)
                   {
                      if(($$seqseq[$i][4] eq "I1R2" && $$seqseq[$j][4] eq "LINK") || ($$seqseq[$j][4] eq "I1R2" && $$seqseq[$i][4] eq "LINK"))
                         {
                           my $iseq="$$seqseq[$i][2]"."$$seqseq[$i][3]";
                           my $jseq="$$seqseq[$j][2]"."$$seqseq[$j][3]";
                           my $distance = distance($iseq,$jseq); 
                           if ($distance == 1)
                              { 
                                 if($$seqseq[$i][4] eq "I1R2") { $possible_error{$$seqseq[$i][5]} = $$seqseq[$j][5]; $if{$$seqseq[$j][5]} = "$$seqseq[$j][2]\t$$seqseq[$j][3]\t$$seqseq[$j][4]\t$$seqseq[$j][5]"; }
                                else{$possible_error{$$seqseq[$j][5]} = $$seqseq[$i][5]; $if{$$seqseq[$i][5]} = "$$seqseq[$i][2]\t$$seqseq[$i][3]\t$$seqseq[$i][4]\t$$seqseq[$i][5]";}
                              }
                         }
                   }
               }
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
