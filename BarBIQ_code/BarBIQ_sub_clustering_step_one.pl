#! /usr/bin/env perl
##########################################################################################################################################
######Description of this code#####
## This code is used to cluster 16S rRNA sequences (I1 and R2) by nucleotide-sequence-clusterizer with distance 3. 
## As an option, you can change the distance by --d No.
## nucleotide-sequence-clusterizer can be download from http://kirill-kryukov.com/study/tools/nucleotide-sequence-clusterizer/.
## please use the output file of BarBIQ_sub_clean_end_primer_and_barcode.pl/BarBIQ_sub_clean_end_primer_and_barcode_single.pl as the inputfile of this code.
##########################################################################################################################################
######how to run this code #####
###command##
## BarBIQ_sub_clustering_step_one.pl --in cluster_I1_R2_file --out outputfile (--d D)
###interpretation##
## --in: the inputfile which is generated from BarBIQ_sub_clean_end_primer_and_barcode.pl/BarBIQ_sub_clean_end_primer_and_barcode_single.pl, so please run that code first.
## --out: outputfile, please set a name for your outputfile
## --d: sequences separated by D or fewer substitutions (the default d value is 3. you can change it using option --d, the largest No. could be 3)
##########################################################################################################################################
#####Install#####
## Please install the nucleotide-sequence-clusterizer and add it to Environment variable (let it possible to be called directly)
## please install the perl Module IPC::System::Simple qw(system) before using this code
##########################################################################################################################################

#####code#####
use strict;
use warnings;
use IPC::System::Simple qw(system); ## for calling external program

#print "Now you are runing $0\n";
#print "The parameters are: @ARGV\n";

##read command##
my ($i,$inputfile,$outputfile,$d,$software,$log_file);
$d=3; # the defult value
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--out") {$outputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--d") {$d = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--soft") {$software = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--log") {$log_file= $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in cluster_I1_R2_file --out outputfile --soft nucleotide-sequence-clusterizer (--log log_file_name --d 3)\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: inputfile\"\n $!";}
if(!$outputfile)  {die "Your input is wrong!!!\n Please input \"--out: outputfile\"\n $!";}
if(!$log_file)  {
     $log_file = $outputfile."_log";
     }
open (my $LOG, '>>', $log_file);
select $LOG;
print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";
print "You are using distance: $d\n";

if(!$software) {die "Your input is wrong!!!\n Please input \"--soft nucleotide-sequence-clusterizer\"\n $!";}
if(-e $outputfile){die "Your output file $outputfile is already existed!!! please check!!!\n $!";}
##read command##

##check the inputfile##
my $single_end="NO";
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
if($info[4] eq "No" && $info[5] eq "No")
   {
    print "This is an single end data!!!\n";
    $single_end="YES";
   }
close FILE;
print "Inputfile is OK!\nStart to calculating:\n";
##check the inputfile##

##sorting the inputfile by cluster ID##
print "sorting...\n";
my $inputfile_sort="$inputfile"."_sort"; ## a file for sorted inputfile by cluster ID
unlink $inputfile_sort; 
my (@out,$out,$string);
open(IN, $inputfile) || die "canot open input file '$inputfile' $!";
@out = sort {$a->[0] cmp $b->[0]}
       map [(split)], <IN>;
open(SORT, '>>', $inputfile_sort) or die "Could not open file '$inputfile_sort' $!";
for $out(@out)
   {
    $string=join("\t",@$out);
    print SORT ("$string\n");
   }
close SORT;
close IN;
print "File sorting finished!!!\n";
##sorting the inputfile by cluster ID##

##Get one cluster's data and do sub-clusterizing##
open (FILE,$inputfile_sort) or die "Could not open file '$inputfile_sort' $!";
open (OUTF,'>>', $outputfile) or die "Could not open file '$outputfile' $!";
chomp($gi=<FILE>);
@info=split(/\s+/,$gi);
my $cluster_name = $info[0];
my @seq_I1=split("", $info[2]); my $I1_seq_length=($#seq_I1+1); ## count the length of I1
my @seq_R2=split("", $info[4]); my $R2_seq_length=($#seq_R2+1); ## count the length of R2
my $t;
if($single_end eq "YES")
    {$t=$I1_seq_length;}
elsif($single_end eq "NO")
    {$t=$I1_seq_length+$R2_seq_length;}
else{die "$single_end is not YES or NO!!! please check!!!\n $!";} ## --t for nucleotide-sequence-clusterizer
my $T="";
for($i=0; $i<$t; $i++)
   {$T="$T".".";} 
my %cluster;
my $outputfile_c = "$outputfile"."_$cluster_name".".fasta";  ##reads of one cluster, for the input of nucleotide-sequence-clusterizer
my $outputfile_o = "$outputfile"."_$cluster_name"."_sub";  ## for the output of nucleotide-sequence-clusterizer
unlink $outputfile_c;
unlink $outputfile_o;
open (OUTFC,'>>', $outputfile_c) or die "Could not open file '$outputfile_c' $!";
if($single_end eq "YES")
   {
    print OUTFC (">$info[1]\n$info[2]\n");
   }
elsif($single_end eq "NO")
   {
    print OUTFC (">$info[1]\n$info[2]$info[4]\n");
   }
else{die "$single_end is not YES or NO!!! please check!!!\n $!";}
my $seqqual="$info[2]\t$info[3]\t$info[4]\t$info[5]";
$cluster{$info[1]}=$seqqual;
while($gi=<FILE>)
   {
    chomp $gi;
    @info=split(/\s+/,$gi);
    if ($info[0] eq $cluster_name) {
       if($single_end eq "YES")
           {
            print OUTFC (">$info[1]\n$info[2]\n");
           }
       elsif($single_end eq "NO")
           {
            print OUTFC (">$info[1]\n$info[2]$info[4]\n");
           }
       else{die "$single_end is not YES or NO!!! please check!!!\n $!";}
       $seqqual="$info[2]\t$info[3]\t$info[4]\t$info[5]";
       $cluster{$info[1]}=$seqqual;
      } else {
       close OUTFC;
       &subcluserizing($outputfile_c,$outputfile_o,$cluster_name);
       unlink $outputfile_c;
       unlink $outputfile_o;
       $cluster_name = $info[0];
       $outputfile_c = "$outputfile"."_$cluster_name".".fasta";  ##reads of one cluster, for the input of nucleotide-sequence-clusterizer
       $outputfile_o = "$outputfile"."_$cluster_name"."_sub";  ## for the output of nucleotide-sequence-clusterizer
       unlink $outputfile_c;
       unlink $outputfile_o;
       undef %cluster;
       open (OUTFC,'>>', $outputfile_c) or die "Could not open file '$outputfile_c' $!";
       if($single_end eq "YES")
           {
            print OUTFC (">$info[1]\n$info[2]\n");
           }
       elsif($single_end eq "NO")
           {
            print OUTFC (">$info[1]\n$info[2]$info[4]\n");
           }
       else{die "$single_end is not YES or NO!!! please check!!!\n $!";}
       $seqqual="$info[2]\t$info[3]\t$info[4]\t$info[5]";
       $cluster{$info[1]}=$seqqual;
      }
    }
       close OUTFC;
       &subcluserizing($outputfile_c,$outputfile_o,$cluster_name);
       unlink $outputfile_c;
       unlink $outputfile_o;
close FILE;


#nucleotide-sequence-clusterizer#
sub subcluserizing 
    {
     system $software, '-i', $_[0], '-o', $_[1], '-d', $d, '-t', $T; ## sub-clustering using nucleotide-sequence-clusterizer
     open (FILENO,$_[1]) or die "Could not open file '$_[1]' $!";
     my $gi=<FILENO>;
     my $sub_cluster;
     if ($gi =~ m{\A>([0-9]*)}) {$sub_cluster=$1;} else {die "Outputfile $_[1] of nucleotide-sequence-clusterizer is not correct!!! $!";}
     while($gi=<FILENO>)
         {
          chomp $gi;
          if ($gi =~ m{\A>([0-9]*)}) 
              {
                $sub_cluster=$1;
              } 
          else{
               my @line=split(/\s+/,$gi);
               print OUTF ("$_[2]\t$sub_cluster\t$line[0]\t$cluster{$line[0]}\n");
              }
         }
     close FILENO;
    }
#nucleotide-sequence-clusterizer#
close OUTF;
unlink $inputfile_sort;
print "Done!!!\n";
##Get one cluster's data and do sub-clusterizing##

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.002
#2022.10.14
