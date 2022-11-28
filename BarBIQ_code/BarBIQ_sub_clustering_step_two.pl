#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used to cluster 16S rRNA sequences (I1 and R2) based on a single position of the reads; and also to generate a representitive sequence (RepSeq) for each sub-cluster(SCluster).
## This step is after the BarBIQ_sub_clustering_step_one.pl.
## and use the output file from BarBIQ_sub_clustering_step_one.pl as the inputfile of this code.
## Please note that this code only accept the Phred+33 quality system.
## if in one subcluster, there is one postion, all reads's quality is lower than threshold quality, then this subcluster will be deleted.
## if in one subcluster, there is only one read,then this subcluster will be deleted.
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_clustering_step_two.pl --in cluster_I1_R2_file_sub1 --out outputfile (--ratio No. --qual No.)
###explaination##
## --in: the inputfile which is generated from BarBIQ_sub_clustering_step_one.pl.
## --out: outputfile, please set a name for your outputfile
## --ratio: a threshold for clustering (the default value is 0.75. you can change it using option --threshold, the value should be beteen 0 to 1)
## --qual: a threshold for considering the quality of each base, when the sequencing quality score was < threshold, the converted score was 0, and when the sequencing quality score was >= threshold, the converted score was equal to the sequencing quality score divided by 41 (highest sequencing quality score) (the default value is 15. you can change it using option --qual, the value should be beteen 0 to 41)
##########################################################################################################################################
#####Install#####
## please install the perl Modules IPC::System::Simple qw(system) and List::Util qw(max) before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system);
use List::Util qw(max);

#print "Now you are runing BarBIQ_sub_clustering_step_two.pl\n";
#print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$outputfile,$threshold_ratio,$threshold_quality, $log_file);
$threshold_ratio=0.75; # the defult value
$threshold_quality=15; # the defult value
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--ratio") {$threshold_ratio = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--qual") {$threshold_quality = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in cluster_I1_R2_file_sub1 --out outputfile (--log log_file_name --ratio 0.75 --qual 15)\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out: outputfile\"\n $!";}
if(!$log_file)
     {
        $log_file = $outputfile."_log";
     }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
print "Your are using the threshold_ratio of $threshold_ratio and the threshold_quality of $threshold_quality\n";

if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
##read command##

##check the inputfile##
if(!(-e $inputfile)) {die "Your input file $inputfile is not existed!!! please check!!!\n $!";}
open (FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
print "Your inputfiles is:\n$inputfile\n";
my $gi=<FILE>;
chomp $gi;
my @info=split(/\s+/,$gi);
my @I1_read=split("", $info[3]);
my $I1_length=($#I1_read+1); ## how many bases of I1
if(!(($#info == 6) && ($info[0] =~ /\Acluster_/)))
   {
    die "Your input file $inputfile is wrong!!! please check!!!\n $!";
   }
my $single_end="NO";
if($info[5] eq "No" && $info[6] eq "No")
    { 
      print "This this a single end data!!!\n";
      $single_end="YES";
    }
close FILE;
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

###Main code###
my $middlefile1="$outputfile"."_s"; ## a middle file for store the subclusters in each round
unlink $middlefile1;
my $middlefile2 = "$middlefile1"."_b"; ## a middle file for store the unfinal subclusters in each round
unlink $middlefile2;
system "cp $inputfile $middlefile1";
while(-e $middlefile1)
    {
     &sub_clusterizing_step_two($middlefile1, $middlefile2);
     unlink $middlefile1;
     if(-e $middlefile2)
         {
              &sort_by_subcluster($middlefile2, $middlefile1);
              unlink $middlefile2;
         }
    }
     
##Get one cluster's data and do sub-clusterizing step two##
sub sub_clusterizing_step_two 
    {
      open (FILESUB2,$_[0]) or die "Could not open file '$_[0]' $!";
      chomp($gi=<FILESUB2>);
      @info=split(/\s+/,$gi);
      my $subcluster_name = "$info[0]"."_$info[1]";
      my $middlefile_c = "$_[0]"."_$subcluster_name";  ##reads of one subcluster
      unlink $middlefile_c;
      open (MIDC,'>>', $middlefile_c) or die "Could not open file '$middlefile_c' $!";
      print MIDC ("$gi\n");
      while($gi=<FILESUB2>)
           {
            chomp $gi;
            @info=split(/\s+/,$gi);
            my $subcluster_name_next="$info[0]"."_$info[1]";
            if ($subcluster_name eq $subcluster_name_next) 
               {
                print MIDC ("$gi\n");
               } else {
                close MIDC;
                &subcluserizing($middlefile_c, $_[1], $subcluster_name);
                unlink $middlefile_c;
                $subcluster_name = $subcluster_name_next;
                $middlefile_c = "$_[0]"."_$subcluster_name";  ##reads of one subcluster
                unlink $middlefile_c;
                open (MIDC,'>>', $middlefile_c) or die "Could not open file '$middlefile_c' $!";
                print MIDC ("$gi\n");  
               }
            }
       close MIDC;
       &subcluserizing($middlefile_c, $_[1], $subcluster_name);
       unlink $middlefile_c;
       close FILESUB2;
#subclustering#
sub subcluserizing 
    {
     open(FILESUB, $_[0]) or die "cannot open input file '$_[0]' $!";
     my ($gi, @info, $seq, @seq, @sequence, %qual, %subcluster1, %subcluster2);
     while($gi=<FILESUB>)
        {
         chomp $gi;
         @info=split(/\s+/,$gi);
         if ($single_end eq "NO")
             {
              $seq="$info[3]"."$info[5]";
              $qual{$info[2]}="$info[4]"."$info[6]";
             }
         elsif ($single_end eq "YES")
             {
              $seq="$info[3]";
              $qual{$info[2]}="$info[4]";
             }
         else {die "$single_end is not YES or NO!!! please check!!!\n $!";}
         @seq=split("",$seq);
         unshift @seq, $info[2];
         push @sequence, [@seq];
         $subcluster1{$info[2]}="$info[0]\t$info[1]";
         $subcluster2{$info[2]}="$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]";
        }
     close FILESUB;
     my ($i, $j, @quality, $q, @A, @T, @C, @G, $qualscore);
     $i=0; # the first read
     @quality=split("",$qual{$sequence[$i][0]});# quality of the read
     $q=0; # the first base;
     for ($j=1;$j<=$#{$sequence[$i]}; $j++)
         {
          if ($sequence[$i][$j] eq "A") { $A[$j] = (ord($quality[$q])-33)/41; if($A[$j]<($threshold_quality/41)) {$A[$j]=0;} $T[$j]=0;$C[$j]=0;$G[$j]=0;$q++;}
       elsif ($sequence[$i][$j] eq "T") { $T[$j] = (ord($quality[$q])-33)/41; if($T[$j]<($threshold_quality/41)) {$T[$j]=0;} $A[$j]=0;$C[$j]=0;$G[$j]=0;$q++;}
       elsif ($sequence[$i][$j] eq "C") { $C[$j] = (ord($quality[$q])-33)/41; if($C[$j]<($threshold_quality/41)) {$C[$j]=0;} $A[$j]=0;$T[$j]=0;$G[$j]=0;$q++;}
       elsif ($sequence[$i][$j] eq "G") { $G[$j] = (ord($quality[$q])-33)/41; if($G[$j]<($threshold_quality/41)) {$G[$j]=0;} $A[$j]=0;$T[$j]=0;$C[$j]=0;$q++;}
       elsif ($sequence[$i][$j] eq "N") { $A[$j]=0;$T[$j]=0;$C[$j]=0;$G[$j]=0;$q++;}
       else                             { print "Your data has this sequence: $sequence[$i][$j], please check!!!\n"; last;}
         }
     for ($i=1; $i<=$#sequence; $i++)
         {
          @quality=split("",$qual{$sequence[$i][0]}); # quality of the read
          $q=0; # the first base;
          for ($j=1;$j<=$#{$sequence[$i]}; $j++)
             {
              if ($sequence[$i][$j] eq "A") { $qualscore = (ord($quality[$q])-33)/41; if($qualscore<($threshold_quality/41)) {$qualscore=0;} $A[$j] = $A[$j]+$qualscore;$q++;}
           elsif ($sequence[$i][$j] eq "T") { $qualscore = (ord($quality[$q])-33)/41; if($qualscore<($threshold_quality/41)) {$qualscore=0;} $T[$j] = $T[$j]+$qualscore;$q++;}
           elsif ($sequence[$i][$j] eq "C") { $qualscore = (ord($quality[$q])-33)/41; if($qualscore<($threshold_quality/41)) {$qualscore=0;} $C[$j] = $C[$j]+$qualscore;$q++;}
           elsif ($sequence[$i][$j] eq "G") { $qualscore = (ord($quality[$q])-33)/41; if($qualscore<($threshold_quality/41)) {$qualscore=0;} $G[$j] = $G[$j]+$qualscore;$q++;}
           elsif ($sequence[$i][$j] eq "N") { $q++;}
           else                    { print "Your data has this sequence: $sequence[$i][$j], please check!!!\n"; last;}
             }
         }
     my ($A, $T, $C, $G, $r, $qualdel, @quality_max, @seq_max, $max_score_2, $max_score, $sum, $doublepos, $max_score_max,$seq_max_max, @score_relative, $doublepos_max, @quality_max_all, $dontsave);
     $doublepos=0;
     $max_score_max=0;
     $dontsave=0;
     for ($r=1; $r<=$#A; $r++)
        {
         $sum = $A[$r]+$T[$r]+$C[$r]+$G[$r];
         if($sum != 0)
             {
              $A=$A[$r]/$sum;
              $T=$T[$r]/$sum;
              $C=$C[$r]/$sum;
              $G=$G[$r]/$sum;
              $max_score = max ($A, $T, $C, $G);
              if ($max_score == $A) {$seq_max[$r] = "A"; @score_relative = ($T, $C, $G);}
           elsif ($max_score == $T) {$seq_max[$r] = "T"; @score_relative = ($A, $C, $G);}
           elsif ($max_score == $C) {$seq_max[$r] = "C"; @score_relative = ($A, $T, $G);}
           elsif ($max_score == $G) {$seq_max[$r] = "G"; @score_relative = ($A, $T, $C);}
           else                     {die; print "something is wrong!!\n"}
              $max_score_2 = max @score_relative;
              $quality_max[$r] = $max_score_2/$max_score;
              if ($quality_max[$r] > $threshold_ratio) 
                 {
                  $doublepos=$r;
                  if ($max_score_max < $quality_max[$r]) {$max_score_max=$quality_max[$r]; $seq_max_max=$seq_max[$r]; $doublepos_max=$doublepos;}
                 }
            } else {$dontsave = 1;last;}
         }
       my (@seq_max_I1, @seq_max_R2,$seq_max_I1, $seq_max_R2);
       if ($dontsave == 0)
          {
           if (!$max_score_max) 
               {
                @seq_max_I1=@seq_max;
                shift @seq_max_I1;
                @seq_max_R2=splice @seq_max_I1, $I1_length;
                $seq_max_I1 = join ("", @seq_max_I1);
                if(@seq_max_R2) 
                    {
                     $seq_max_R2 = join ("", @seq_max_R2);
                    }
                else{
                     $seq_max_R2 = "No";
                    }
                my $reads=$#sequence+1;
                if($reads>1)
                   {
                    open(OUTF, ">>", $outputfile) or die "cannot open input file '$outputfile' $!";
                    print OUTF ("$_[2]\t$reads\t$seq_max_I1\t$seq_max_R2\n");
                    close OUTF;
                   }
               }
           else{ 
                open(OUTF2, ">>", $_[1]) or die "cannot open input file '$_[1]' $!";
                for ($i=0;$i<=$#sequence; $i++)
                     {
                      if($sequence[$i][$doublepos_max] eq $seq_max_max)
                          {
                           print OUTF2 ("$subcluster1{$sequence[$i][0]}\_1\t$subcluster2{$sequence[$i][0]}\n");
                          }
                      else{
                           print OUTF2 ("$subcluster1{$sequence[$i][0]}\_2\t$subcluster2{$sequence[$i][0]}\n");
                          }
                      }
                close OUTF2;
               }
           }
    }
#subclustering#
}
#sort by subclusters#
sub sort_by_subcluster
    {
     my (@out, $out, $string);
     open(IN, $_[0]) || die "canot open input file '$_[0]' $!";
     @out = sort {$a->[0] cmp $b->[0]
               or $a->[1] cmp $b->[1]}
            map [(split)], <IN>;
     open(SORT, '>>', $_[1]) or die "Could not open file '$_[1]' $!";
     for $out(@out)
        {
         $string=join("\t",@$out);
         print SORT ("$string\n");
        }
     close SORT;
     close IN;
    }
#sort by subclusters#
print "Done!!!\n";
##Finshed Get one cluster's data and do sub-clusterizing##
###Main code###
####end####

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
