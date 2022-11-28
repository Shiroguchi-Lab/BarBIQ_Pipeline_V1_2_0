#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
##This code is used to remove contaminated bacteria using a control.
##An information file to describe the samples and controls shoud be prepared, see example IF.txt
##########################################################################################################################################
######how to run this code #####
###command##
#BarBIQ_final_review_contamination_bac_123rep.pl --in file1 --IF file2 (--Factor No.)
###explaination##
#file1: Out put file from BarBIQ_final_compare_bacteria_count.pl
#file2: An information file to describe the samples and controls, see example IF.txt
#Factor No.: is an option, if you used different number of droplets to prepare library for the sample and control, you may need to use this factor to normalize.
##########################################################################################################################################
#####Install#####
##please install the perl Module Math::Round before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
# use Statistics::Basic qw(:all);
use Math::Round;

# print "Now you are runing $0\n";
# print "The parameters are: @ARGV\n";


##read command##
my ($i,$inputfile,$IF_file, $log_file);
my $factor = 1;
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--IF") {$IF_file = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--Factor") {$factor = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in inputfile --IF IF_file (option: --Factor 2.5, 1 is default --log log_file)\"\n $!";}
    }
if(!$log_file) {
                 $log_file = $inputfile."_log";
               }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
my $outputfile = "$inputfile"."_ave_del_contamination";
my $outputfile_2 = "$inputfile"."_raw_del_contamination";
unlink $outputfile;
unlink $outputfile_2;
if(!$IF_file)  {die "Your input is wrong!!!\n Please input \"--IF file\"\n $!";}
print "You are using size factor $factor for the Blank control data, if the factor is not 1 means your used different number of droplets for the sample measurement and blank!!!\n";
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfiles is:\n$inputfile\n";
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
if(!($info[0] eq "ID"))
   {
    die "Your input file $inputfile is wrong!!! please check!!!\n $!";
   }
close FILE;
print "Inputfiles are OK!\nStart to calculating:\n";
##check the inputfile##

###Main code### 
open(IF, $IF_file) or die "canot open input file '$IF_file' $!";
$gi=<IF>;
chomp $gi;
@info=split(/\s+/,$gi);
my $index;
my $sample;
my $type;
for (my $i=0; $i<=$#info; $i++)
    {
        if($info[$i] eq "Index") {$index=$i;}
     elsif($info[$i] eq "Sample") {$sample=$i;}
     elsif($info[$i] eq "Type") {$type=$i;}
     else {die "Your input file $IF_file is wrong!!! please check!!!005\n $!";}
    }
my %types;
my $sample_no=0;
my @Blank;
while($gi=<IF>)
    {
      chomp $gi;
      @info=split(/\s+/,$gi);
      if(exists $types{$info[$type]})
         {
           $types{$info[$type]} = $types{$info[$type]}."\t$info[$index]";
         }
     else{ $types{$info[$type]} = $info[$index];}
     if($info[$type] eq "Blank")
        {
          push @Blank, $info[$index];
        }
     $sample_no++;
    }
close IF;

open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
$gi=<FILE>;
chomp $gi;
@info=split(/\s+/,$gi);
my %indexs;
my %sums;
for (my $i=1; $i<=$sample_no; $i++)
       {
         $indexs{$info[$i]} = $i;
         $sums{$info[$i]}=0;
       }
my @data;

while($gi=<FILE>)
    {
       chomp $gi;@info=split(/\s+/,$gi);
       push @data, [@info];
       foreach my $key (keys %indexs)
         {
           $sums{$key}=$sums{$key}+$info[$indexs{$key}];
         }
    }
close FILE;

## Normalization start ##
foreach my $key (sort keys %types)
    {
      if ($key ne "Blank")
          {
      my @indexindex=split(/\s+/, $types{$key});
      my $max=$sums{$indexindex[0]};
      if($#indexindex) 
         {
           for(my $in=1; $in<=$#indexindex; $in++)
              {
                if($sums{$indexindex[$in]}>$max) { $max=$sums{$indexindex[$in]}; }
              }
         }
      &Nomalization($max,\@indexindex)
        }
      else {
             my @indexindex=split(/\s+/, $types{$key});
             &Nomalization_factor($factor,\@indexindex) 
           }
    }

sub Nomalization_factor
   {
     my ($f, $indexindexindex) = @_;
     for(my $i=0; $i<=$#data; $i++)
         {
            for(my $j=0; $j<=$#$indexindexindex; $j++)
                {
                  $data[$i][$indexs{$$indexindexindex[$j]}]=$data[$i][$indexs{$$indexindexindex[$j]}]*$f;
                }
         }
   }

sub Nomalization
   {
     my ($maxmax, $indexindexindex) = @_;
     for(my $i=0; $i<=$#data; $i++)
         {
            for(my $j=0; $j<=$#$indexindexindex; $j++)
                {
                  $data[$i][$indexs{$$indexindexindex[$j]}]=$data[$i][$indexs{$$indexindexindex[$j]}]*$maxmax/$sums{$$indexindexindex[$j]};
                }
         }
   }

## Normalization done##

## Delete the contamination start ##

my @data_del_con;
my $ppp=1;
foreach my $key (sort keys %types)
    {
      my @indexindex=split(/\s+/, $types{$key});
      if ($key ne "Blank")
         {
           print "Sample name: $key\n";
           $data_del_con[0][$ppp]=$key;
           &Del_Con(\@indexindex, $key);
           $ppp++;
         }
    }
sub Del_Con
     {
       my ($indexindexindex, $typetype) = @_;
       my @data_del_con_raw;
       $data_del_con_raw[0][0]="ID";
       for(my $j=0; $j<=$#$indexindexindex; $j++)
                {
                 $data_del_con_raw[0][$j+1]=$$indexindexindex[$j];
                }
       for(my $i=0; $i<=$#data; $i++)
         {
            my @count;
            my @Blank_count;
            $data_del_con_raw[$i+1][0]=$data[$i][0];
            for(my $j=0; $j<=$#$indexindexindex; $j++)
                {
                  push @count, $data[$i][$indexs{$$indexindexindex[$j]}];
                }
            for(my $j=0; $j<=$#Blank; $j++)
               {
                   push @Blank_count, $data[$i][$indexs{$Blank[$j]}];
               }
      my $mean_s;
      my $sd_s;
      my $mean_b;
      my $sd_b;
        if ($#count)
           {
            my $rps = $#count +1;
           if($i == 0) {
            print "You only have $rps repeat of this sample, so using measued sd!\n";
             }
            $mean_s=&mean(\@count); $mean_s=~s/,//i;
            $sd_s=&stddev(\@count); $sd_s=~s/,//i;
           }
        elsif ($#count == 0)
           {
            if($i == 0) {
            print "You only have one repeat of this sample, so using sampling nosie as sd!\n";}
            $mean_s=$count[0];$mean_s=~s/,//i;
            $sd_s=sqrt($count[0]); $sd_s=~s/,//i;
           }
        if ($#Blank_count >1)
          {
             my $rps = $#Blank_count +1;
            if($i == 0) {
            print "You only have $rps repeats of the blank, so using measued sd!\n";}
            $mean_b=&mean(\@Blank_count); $mean_b=~s/,//i;
            $sd_b=&stddev(\@Blank_count); $sd_b=~s/,//i;
             $sd_b=$sd_b*3.27;
            # if ( $sd_b < $mean_b * 0.1 ) { $sd_b = $mean_b * 0.1 }
           }
       elsif ($#Blank_count == 0)
           {
            if($i == 0) {
            print "You only have one repeat of the blank, so using sampling nosie as sd!\n";}
            $mean_b=$Blank_count[0]; $mean_b=~s/,//i;
            $sd_b=$mean_b * 0.3; $sd_b=~s/,//i;
           }
       elsif ($#Blank_count == 1)
           {
             if($i == 0) {
            print "You only have 2 repeats of the blank, so using measued sd!\n";}
            $mean_b=&mean(\@Blank_count); $mean_b=~s/,//i;
            $sd_b=&stddev(\@Blank_count); $sd_b=~s/,//i;
            $sd_b=$sd_b*3.27;
            if ($sd_b < ($mean_b * 0.1)) { $sd_b = $mean_b * 0.1 }
           }
          if($mean_b == 0 )
             {
                $data_del_con[$i+1][$ppp] = $mean_s;
                for(my $r=0; $r<=$#count; $r++)
                    {
                     $data_del_con_raw[$i+1][$r+1] = $count[$r];
                    }
             }
         else{
            if (($mean_s-$sd_s) > ($mean_b+$sd_b))
               {
                 $data_del_con[$i+1][$ppp] = round($mean_s-$mean_b);
                 for(my $r=0; $r<=$#count; $r++)
                    {
                     $data_del_con_raw[$i+1][$r+1] = round($count[$r]-$mean_b);
                    }
               }
           else{
                 $data_del_con[$i+1][$ppp] = 0;
                 for(my $r=0; $r<=$#count; $r++)
                    {
                     $data_del_con_raw[$i+1][$r+1] = 0;
                    }
               }
         }
        }
  
      my $outputfile_3 = "$outputfile_2"."_$typetype";
      unlink $outputfile_3;
      open(OUTF3, '>>', $outputfile_3) or die "canot open output file '$outputfile_3' $!";
      for my $out(@data_del_con_raw)
           {
            my $string=join("\t",@$out);
            print OUTF3 ("$string\n");
           }
      close OUTF3;
#### You can output the raw count for each replicate by removing the following line.
     unlink $outputfile_3;
     }
## Delete the contamination done## 
for(my $i=0; $i<=$#data; $i++)
   {
     $data_del_con[$i+1][0]=$data[$i][0];
   }
$data_del_con[0][0]="ID";
## save data ##
open(OUTF, '>>', $outputfile) or die "canot open output file '$outputfile' $!";
for my $out(@data_del_con)
  {
  my $string=join("\t",@$out);
  print OUTF ("$string\n");
  }
close OUTF;
print "Done!!!\n";

sub mean {
    my $data = $_[0];
    my ($sum)=0;
    my ($sqsum)=0;
    foreach my $x(@{$data}) {
       $sum += $x;
       }
    my $mean=$sum/($#{$data} +1);
    return $mean;
    }

sub stddev {
    my $data = $_[0];
    if(@$data == 1){
        return 0;
     }
    my $mean = &mean($data);
    my $sqtotal = 0;
    foreach(@$data) {
          $sqtotal += ($mean-$_) ** 2;
        }
   my $std = ($sqtotal / (@$data-1)) ** 0.5;
   return $std;
  }
       
###Main code###
####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.20
