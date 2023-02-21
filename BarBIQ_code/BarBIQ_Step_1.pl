#! /usr/bin/env perl
##########################################################################################################
#####Description of this code#####
# This code is the first step of the BarBIQ pipline.
# It is a combined pipeline of  
#    BarBIQ_sub_clean_quality_R1.pl
#    BarBIQ_sub_barcode_fix.pl
#    BarBIQ_sub_barcode_clustering.pl
#    BarBIQ_sub_index_and_leakage.pl
#    BarBIQ_sub_clean_I1R2_by_R1.pl
#    BarBIQ_sub_average_quality_of_each_position.pl 
#    BarBIQ_sub_statistic_reads_per_barcode_by_R1.pl
#    BarBIQ_M_I1R2.pl 
# So please make sure these codes are all under the same directory of this combined code.
# If you use this combined code instead of those separated codes respectively, meaning you want to use all default parameters without any modification.
# The default parameters please check those codes in detail respectively. 
# The input data should be fastq format (.fastq or .fastq.gz) files from Miseq directly, and all samples measured in the same MiSeq run should be analyzed together. 
# Output files are separated for each index.
# This code uses these input files(see examples), please prepare in advance:
# BarBIQ_example_inputfile.txt : the data file names
# BarBIQ_example_parameter.txt : necessary parameters
###########################################################################################################
#####how to run this code #####
##command##
## BarBIQ_Step_1.pl --in BarBIQ_example_inputfile.txt --p BarBIQ_example_parameter.txt (--middle Yes)
##interpretation##
# --in: a file should be prepared like the example BarBIQ_example_inputfile.txt 
# --p: a file including necessary parameters should be prepared like the example BarBIQ_example_parameter.txt
# --middle: if you want to keep the middlefiles which generated during processing, please set this parameter as "Yes", default is "No"
########################################################################################################## 
#####Install#####
## Please install the perl modules: Bio::SeqIO, Bio::Seq, File::Copy, File::Path qw(make_path), and IPC::System::Simple qw(system) before use this code
## Please install the nucleotide-sequence-clusterizer and add it to Environment variable(let it possibale to be called directly)
############################################################################################################

use strict;
use warnings;
use IPC::System::Simple qw(system);
use File::Path qw(make_path);
use File::Copy;
use Parallel::ForkManager;

# print "Now you are running program: $0\n";
# print "Your parameters are: @ARGV\n";
# print "Started at: ";
# print scalar localtime;
# print "\n";

##read command##
my ($i,$inputfile,$parefile);
my $keep_middle="No";
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--p") {$parefile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--middle") {$keep_middle = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in BarBIQ_example_inputfile.txt --p BarBIQ_example_parameter.txt --out outputfile\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: BarBIQ_example_inputfile.txt\"\n $!";}
if(!$parefile)  {die "Your input is wrong!!!\n Please input \"--p BarBIQ_example_parameter.txt\"\n $!";}
##read command##

##Check the parameter file##
my ($I1_length, $R2_length, $MAX_PROCESSES, $nucleotide);
my (@real_data, %index, @std_data, $fixed_barcode_file, %samplename, %R2files, %I1files, %R1files, %Sampletypes, $fixed_barcode_std_file, $gi, @info, $pippath, $datapath, $code);
open(FILE,$parefile) or die "Could not open file '$parefile' $!"; # open inputfile
while(my $gi=<FILE>)  ## read Real-data file names
    {
     chomp $gi;
     my @info=split(/\s+/,$gi);
     if($info[0] eq "Sample_fix_base1:" || $info[0] eq "Sample_fix_base2:" || $info[0] eq "Sample_fix_base3:" || $info[0] eq "Sample_fix_base4:")
          {
          }
  elsif($info[0] eq "SpikeIn_fix_base1:" || $info[0] eq "SpikeIn_fix_base2:" || $info[0] eq "SpikeIn_fix_base3:" || $info[0] eq "SpikeIn_fix_base4:")
          {
          }
  elsif($info[0] eq "Perl_code_path:") {$pippath = $info[1];}
  elsif($info[0] eq "Data_path:") {$datapath = $info[1];}
  elsif($info[0] eq "I1_length:") {$I1_length = $info[1];}
  elsif($info[0] eq "R2_length:") {$R2_length = $info[1];}
  elsif($info[0] eq "Max_processes:") {$MAX_PROCESSES = $info[1];}
  elsif($info[0] eq "nucleotide-sequence-clusterizer:") {$nucleotide = $info[1];}
  elsif($info[0] ne "#")
      {
        print ("$info[0]\n");
      }
  }
close FILE;
my @parts = split('/', $datapath);
my $outputfile1=$parts[$#parts];
my $log_file = "$datapath"."$outputfile1"."_log";
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are running program: $0\n";
print "Your parameters are: @ARGV\n";
print "Started at: ";
print scalar localtime;
print "\n";

open(FILE,$parefile) or die "Could not open file '$parefile' $!"; # open inputfile
$fixed_barcode_std_file = "$datapath"."BarBIQ_example_fixed_base_std.txt";
unlink $fixed_barcode_std_file;
open (FIXSTD,'>>',$fixed_barcode_std_file) or die "Could not open file '$fixed_barcode_std_file' $!"; # open file
$fixed_barcode_file = "$datapath"."BarBIQ_example_fixed_base.txt";
my $fixsample = 0;
my $fixstd = 0;
unlink $fixed_barcode_file;
open (FIX,'>>',$fixed_barcode_file) or die "Could not open file '$fixed_barcode_file' $!"; # open file
while($gi=<FILE>)  ## read Real-data file names
    {
     chomp $gi;
     @info=split(/\s+/,$gi);
     if($info[0] eq "Sample_fix_base1:" || $info[0] eq "Sample_fix_base2:" || $info[0] eq "Sample_fix_base3:" || $info[0] eq "Sample_fix_base4:")
          {
             print FIX ("$info[1]\n");
             $fixsample++;
          }
  elsif($info[0] eq "SpikeIn_fix_base1:" || $info[0] eq "SpikeIn_fix_base2:" || $info[0] eq "SpikeIn_fix_base3:" || $info[0] eq "SpikeIn_fix_base4:")
          {
             print FIXSTD ("$info[1]\n");
             $fixstd++;
          }
  }
close FILE;
close FIX;
close FIXSTD;
if(!$pippath)  {die "Your parameter file is wrong!!\n Please privide Perl_code_path\n $!";}
if(!$datapath) {die "Your parameter file is wrong!!\n Please privide Data_path\n $!";}
if(!$MAX_PROCESSES) {die "Your parameter file is wrong!!\n Please privide MAX_PROCESSES\n $!";}
if($fixsample == 4 && $fixstd == 4){} else {die "Your parameter file is wrong!!\n Please privide four kinds of fixed bases!!!\n $!";}
my $outputfile = "$datapath"."$outputfile1";
##Finish check the parameter file##


##check the inputfile##
# my (@real_data, %index, @std_data, $fixed_barcode_file, %samplename, %R2files, %I1files, %R1files, %Sampletypes, $fixed_barcode_std_file, $gi, @info, $pippath, $datapath, $code);
open(FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
chomp ($gi=<FILE>);
@info=split(/\s+/,$gi);
my ($indexcol, $R1col, $R2col, $I1col, $Typecol);
for($i=0; $i<=$#info; $i++)
     {  
        if($info[$i] eq "Index"){$indexcol = $i;}
     elsif($info[$i] eq "I1_file"){$I1col = $i;}
     elsif($info[$i] eq "R1_file"){$R1col = $i;}
     elsif($info[$i] eq "R2_file"){$R2col = $i;}
     elsif($info[$i] eq "SampleType"){$Typecol = $i;}
     #else{die "Your input is wrong!!!\n The BarBIQ_example_inputfile.txt should have these columns:Index	I1_file	R1_file	R2_file	SampleType \n $!";}
    }
if(defined $indexcol && defined $R1col && defined $R2col && defined $I1col && defined $Typecol) 
      {
      }
else {
         die "Your input is wrong!!!\n The BarBIQ_example_inputfile.txt should have these columns:Index I1_file R1_file R2_file SampleType \n $!";
      }
while($gi=<FILE>)  ## read Real-data file names
    {
     chomp $gi;
     @info=split(/\s+/,$gi); 
     if($info[$Typecol] eq "SpikeIn")
        { 
          push @std_data, $info[$R1col];
        }
    else{
          push @real_data, $info[$R1col];
        }
       $index{$info[$R1col]}=$info[$indexcol];
       $I1files{$info[$indexcol]}=$info[$I1col];
       $R1files{$info[$indexcol]}=$info[$R1col];
       $R2files{$info[$indexcol]}=$info[$R2col];
       $Sampletypes{$info[$indexcol]}=$info[$Typecol];
       if($info[$R1col] =~ m/_R1_/s)
         {
            $samplename{$info[$indexcol]} = $`;
         }
      else{die "Your input is wrong!!!\n In the BarBIQ_example_inputfile.txt R1_file name should have _R1_ \n $!";}
    }
close FILE; # close inputfile
if(!@real_data) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!@std_data) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%index) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%I1files) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%R1files) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%R2files) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%Sampletypes) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
print "Inputfile is OK!\n";
##Finish check the inputfile##


##Run BarBIQ_sub_clean_quality_R1.pl##
my (@real_data_SBC1,@std_data_SBC1, @real_data_path, @std_data_path);
for($i=0; $i<=$#real_data; $i++)
     {
      $real_data_path[$i] = "$datapath"."$real_data[$i]";
      $real_data_SBC1[$i]="$real_data_path[$i]"."_c.fasta";
      unlink $real_data_SBC1[$i];
      $code="$pippath"."BarBIQ_sub_clean_quality_R1.pl";
      system $code,'--in',$real_data_path[$i],'--out',$real_data_SBC1[$i],'--index',$index{$real_data[$i]},"--log", $log_file;
     }
for($i=0; $i<=$#std_data; $i++)
     {
      $std_data_path[$i]="$datapath"."$std_data[$i]";
      $std_data_SBC1[$i]="$std_data_path[$i]"."_c.fasta";
      unlink $std_data_SBC1[$i];
      $code="$pippath"."BarBIQ_sub_clean_quality_R1.pl";
      system $code,'--in',$std_data_path[$i],'--out',$std_data_SBC1[$i],'--index',$index{$std_data[$i]},"--log", $log_file;
     }
##Finish BarBIQ_sub_clean_quality_R1.pl##

##Run BarBIQ_sub_barcode_fix.pl##
my ($inputfile_names_SBC2, $inputfile_names_std_SBC2, $outputfile_name_SBC2, $outputfile_name_std_SBC2);
$inputfile_names_SBC2="$outputfile"."_sample.txt";
unlink $inputfile_names_SBC2;
open (NAME,'>>',$inputfile_names_SBC2) or die "Could not open file '$inputfile_names_SBC2' $!"; # open file
for($i=0; $i<=$#real_data_SBC1; $i++)
     {
       print NAME ("$real_data_SBC1[$i]\n");
     }
close NAME; # close file
$outputfile_name_SBC2="$outputfile"."_sample.fasta";
unlink $outputfile_name_SBC2;
$inputfile_names_std_SBC2="$outputfile"."_std.txt";
unlink $inputfile_names_std_SBC2;
open (NAME,'>>',$inputfile_names_std_SBC2) or die "Could not open file '$inputfile_names_std_SBC2' $!"; # open file
for($i=0; $i<=$#std_data_SBC1; $i++)
     {
       print NAME ("$std_data_SBC1[$i]\n");
     }
close NAME;# close file
$outputfile_name_std_SBC2="$outputfile"."_std.fasta";
unlink $outputfile_name_std_SBC2;
$code="$pippath"."BarBIQ_sub_barcode_fix.pl";
system $code,'--in',$inputfile_names_SBC2,'--out',$outputfile_name_SBC2,'--fixed',$fixed_barcode_file,"--log", $log_file;
system $code,'--in',$inputfile_names_std_SBC2,'--out',$outputfile_name_std_SBC2,'--fixed',$fixed_barcode_std_file,"--log", $log_file;
##BarBIQ_sub_barcode_fix.pl##

##Run BarBIQ_sub_barcode_clustering.pl##
my ($inputfile_names_SBC3, $outputfile_name_SBC3);
$inputfile_names_SBC3="$outputfile"."_all.txt";
unlink $inputfile_names_SBC3;
open (NAME,'>>',$inputfile_names_SBC3) or die "Could not open file '$inputfile_names_SBC3' $!"; # open file
print NAME ("$outputfile_name_SBC2\n");
print NAME ("$outputfile_name_std_SBC2\n");
close NAME; # close file
$outputfile_name_SBC3="$outputfile"."_clusterizerD2";
unlink $outputfile_name_SBC3;
open (FILE,,$outputfile_name_SBC2) or die "Could not open file '$outputfile_name_SBC2' $!"; # open the inputfile to check the length of read for --t
chomp ($gi=<FILE>);
chomp ($gi=<FILE>);
my @read=split(//,$gi);
my $length=$#read+1;
close FILE; # close file
$code="$pippath"."BarBIQ_sub_barcode_clustering.pl";
system $code,'--in',$inputfile_names_SBC3,'--out',$outputfile_name_SBC3,'--soft',$nucleotide,'--d','2','--t',$length,"--log", $log_file;
##Finish BarBIQ_sub_barcode_clustering.pl##

##Run BarBIQ_sub_index_and_leakage.pl##
my $index_file_SBC4="$outputfile"."_index.txt";
unlink $index_file_SBC4;
open (INDEX,'>>',$index_file_SBC4) or die "Could not open file '$index_file_SBC4' $!"; # open file
print INDEX ("indexes\t");
for($i=0; $i<=$#real_data; $i++)
    {
     print INDEX ("$index{$real_data[$i]}\t");
    }
for($i=0; $i<=$#std_data; $i++)
    {
     print INDEX ("$index{$std_data[$i]}\t"); 
    }
print INDEX ("\n");
close INDEX; # close file;
$code="$pippath"."BarBIQ_sub_index_and_leakage.pl";
system $code,'--in',$outputfile_name_SBC3,'--out',$outputfile,'--index',$index_file_SBC4,"--log", $log_file;
##Finish BarBIQ_sub_index_and_leakage.pl##

##Run BarBIQ_sub_clean_I1R2_by_R1.pl##
my @I1_all;
my @R2_all;
for($i=0; $i<=$#real_data_path; $i++)
    {
     if($real_data_path[$i] =~ m/_R1_/s)
         {
          my $inputfile_name_BarBIQI1=$`."_I1_"."$'";
          my $inputfile_name_BarBIQR2=$`."_R2_"."$'";
          my $inputfile_name_BarBIQcluster="$outputfile"."_$index{$real_data[$i]}";
          if($inputfile_name_BarBIQI1 =~ /.fastq.gz\z/)
             {    
              my $outputfile_name_BarBIQI1_each=$`."_R1.fastq";
              push @I1_all, $outputfile_name_BarBIQI1_each;
             }
       elsif($inputfile_name_BarBIQI1 =~ /.fastq\z/)
             {
              my $outputfile_name_BarBIQI1_each=$`."_R1.fastq";
              push @I1_all, $outputfile_name_BarBIQI1_each;
             }
       else  { die "Something is wrong!005\n";}
          if($inputfile_name_BarBIQR2 =~ /.fastq.gz\z/)
             {  
              my $outputfile_name_BarBIQR2_each=$`."_R1.fastq";
              push @R2_all, $outputfile_name_BarBIQR2_each;
             }
       elsif($inputfile_name_BarBIQR2 =~ /.fastq\z/)
             {
              my $outputfile_name_BarBIQR2_each=$`."_R1.fastq";
              push @R2_all, $outputfile_name_BarBIQR2_each;
             }
               $code="$pippath"."BarBIQ_sub_clean_I1R2_by_R1.pl";
              system $code, '--I1', $inputfile_name_BarBIQI1, '--R2', $inputfile_name_BarBIQR2, '--cluster', $inputfile_name_BarBIQcluster, "--log", $log_file;
         }
    else{ die "Something is wrong!006\n";}
    }
##Finish BarBIQ_sub_clean_I1R2_by_R1.pl##

##Merge all cleaned fastq files##
my $outputfile_name_BarBIQI1="$outputfile"."_I1_all_R1.fastq";
my $outputfile_BarBIQ_I1_ave_qual="$outputfile"."_I1_all_R1_ave_qual";
my $outputfile_BarBIQ_I1_cut="$outputfile"."_I1_cut_position";
my $outputfile_name_BarBIQR2="$outputfile"."_R2_all_R1.fastq";
my $outputfile_BarBIQ_R2_ave_qual="$outputfile"."_R2_all_R1_ave_qual";
my $outputfile_BarBIQ_R2_cut="$outputfile"."_R2_cut_position";
unlink $outputfile_name_BarBIQI1;
unlink $outputfile_name_BarBIQR2;

open my $out, '>>', $outputfile_name_BarBIQI1 or die "Could not open '$outputfile_name_BarBIQI1' for appending\n"; 
foreach my $file (@I1_all) {
    if (open my $in, '<', $file) {
        while (my $line = <$in>) {
            print $out $line;
        }
        close $in;
    } else {
        warn "Could not open '$file' for reading\n";
    }
}
close $out;

open $out, '>>', $outputfile_name_BarBIQR2 or die "Could not open '$outputfile_name_BarBIQR2' for appending\n";
foreach my $file (@R2_all) {
    if (open my $in, '<', $file) {
        while (my $line = <$in>) {
            print $out $line;
        }
        close $in;
    } else {
        warn "Could not open '$file' for reading\n";
    }
}
close $out;
##Finish Merge all cleaned fastq files##


##Run BarBIQ_sub_average_quality_of_each_position.pl##
$code="$pippath"."BarBIQ_sub_average_quality_of_each_position.pl";
system $code, "--in", $outputfile_name_BarBIQI1, "--cut", $outputfile_BarBIQ_I1_cut, "--out", $outputfile_BarBIQ_I1_ave_qual,"--log", $log_file;
system $code, "--in", $outputfile_name_BarBIQR2, "--cut", $outputfile_BarBIQ_R2_cut, "--out", $outputfile_BarBIQ_R2_ave_qual,"--log", $log_file;
##Finish BarBIQ_sub_average_quality_of_each_position.pl##

##Run BarBIQ_add_statistic_reads_per_barcode_by_R1.pl##
open (INXFILE,$index_file_SBC4) or die "Could not open file '$index_file_SBC4' $!"; # open inputfile
my $gigi=<INXFILE>;
chomp($gigi);
my @indexeses=split(/\s+/, $gigi);
if ($indexeses[0] ne "indexes"){die "Your input is wrong!!!\n The index file has not followed the BarBIQ_example_index.txt\n $!";}
for ($i=1; $i<=$#indexeses; $i++) {
   my $R1_clustered_file = "$outputfile"."_$indexeses[$i]";  ##outputfile of $i index
   $code="$pippath"."BarBIQ_sub_statistic_reads_per_barcode_by_R1.pl";
   system $code, "--in", $R1_clustered_file, "--log", $log_file;
   }
close INXFILE;
##Finish BarBIQ_add_statistic_reads_per_barcode_by_R1.pl##

##Delete the middle files ##
if($keep_middle eq "No")
     {
      print "You are deleting the intermediate files\n";
      for($i=0; $i<=$#real_data_SBC1; $i++)
          {
           unlink $real_data_SBC1[$i];
          }
      for($i=0; $i<=$#std_data_SBC1; $i++)
          {
           unlink $std_data_SBC1[$i];
          }
      unlink ($inputfile_names_SBC2,$inputfile_names_std_SBC2, $outputfile_name_SBC2, $outputfile_name_std_SBC2);
      unlink ($inputfile_names_SBC3, $outputfile_name_SBC3);
      unlink $index_file_SBC4;
      unlink ($outputfile_name_BarBIQI1, $outputfile_name_BarBIQR2);
      unlink @I1_all;
      unlink @R2_all;
      unlink ($fixed_barcode_file, $fixed_barcode_std_file);
     }    
##Finish delete the middle files ##

## Move files into each folder for each index##
opendir my $dir, $datapath or die "Cannot open directory: $datapath $!";
my @allfiles = readdir $dir;
closedir $dir;
my $dir_qfn = "$datapath"."/Processed";
make_path($dir_qfn);
# $dir_qfn = "$datapath"."/Processed/SpikeIn";
# make_path($dir_qfn);
foreach my $key (sort keys %Sampletypes)
   {
      print "Moving $key\n";
      if($Sampletypes{$key} eq "SpikeIn")
       {
         if($R1files{$key} =~ m/_R1_/s)
         {
           my @array = grep {/^$`/} @allfiles;
         for($i=0; $i<=$#array; $i++)
            {   
             if(!($array[$i] =~ /.fastq.gz\z/))
                {
                 my $original_file = "$datapath".$array[$i];
                 # my $new_file = "$datapath"."Processed/SpikeIn/".$array[$i];
                 # move($original_file, $new_file);
                 unlink $original_file;
                 # print "Move\t$original_file\tas\t$new_file\n";
                 print "Remove\t$original_file\n";
                }
            }
           my $original_file = "$outputfile"."_$key";
           # my $new_file = "$datapath"."Processed/SpikeIn/".$outputfile1."_$key";
           # move($original_file, $new_file);
           unlink $original_file;
           # print "Move\t$original_file\tas\t$new_file\n";
           print "Remove\t$original_file\n";
          }
        else{die "Something is wrong!015\n";}
       }
     else{
       $dir_qfn = "$datapath"."Processed/".$key;
       make_path($dir_qfn);
       if($R1files{$key} =~ m/_R1_/s)
         {
           my @array = grep {/^$`/} @allfiles;
         for($i=0; $i<=$#array; $i++)
            {
              if(!($array[$i] =~ /.fastq.gz\z/))
                {
                 my $original_file = "$datapath".$array[$i];
                 my $new_file = "$dir_qfn"."/".$array[$i];
                 move($original_file, $new_file);
                 print "Move\t$original_file\tas\t$new_file\n";
                }
            }
           my $original_file = "$outputfile"."_$key";
           my $new_file = "$dir_qfn"."/".$outputfile1."_$key";
           move($original_file, $new_file);  
           print "Move\t$original_file\tas\t$new_file\n";
        }
        else{die "Something is wrong!015\n";}
       }
   }
opendir my $dir2, $datapath or die "Cannot open directory: $datapath $!";
@allfiles = readdir $dir2;
closedir $dir2;
my @array = grep {/^$outputfile1/} @allfiles;
for($i=0; $i<=$#array; $i++)
     {
                 my $original_file = "$datapath".$array[$i];
                 my $new_file = "$datapath"."Processed/".$array[$i];
                 move($original_file, $new_file);
                 print "Move\t$original_file\tas\t$new_file\n";
     }
##Finish move files into each folder for each index##

##Run BarBIQ_M_I1R2.pl##
my ($I1_fastq, $R2_fastq, $outputfile2,$log_file2,$I1_end,$R2_end);
my $outputfile_BarBIQ_I1_cut2=$datapath."Processed/".$outputfile1."_I1_cut_position";
open (FILE,$outputfile_BarBIQ_I1_cut2) or die "Could not open file '$outputfile_BarBIQ_I1_cut2' $!"; # open the inputfile to check the length of read for --t
chomp ($I1_end=<FILE>);
close FILE; # close file
print "You are using the threshold of $I1_end to cut the reads I1\n";
my $outputfile_BarBIQ_R2_cut2=$datapath."Processed/".$outputfile1."_R2_cut_position";
open (FILE,$outputfile_BarBIQ_R2_cut2) or die "Could not open file '$outputfile_BarBIQ_R2_cut2' $!"; # open the inputfile to check the length of read for --t
chomp ($R2_end=<FILE>);
close FILE; # close file
print "You are using the threshold of $R2_end to cut the reads R2\n";
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);
$code = "$pippath"."BarBIQ_M_I1R2.pl";
foreach my $key (sort keys %Sampletypes)
   {
    if($Sampletypes{$key} ne "SpikeIn")
       {
         my $bar_file = $datapath."Processed/".$key."/".$outputfile1."_$key";
         if($I1files{$key} =~ /.fastq.gz\z/)
            {
              $I1_fastq = $datapath."Processed/".$key."/".$`.".fastq";
            }
         else{die "Something is wrong!016\n";}
          if($R2files{$key} =~ /.fastq.gz\z/)
            {
              $R2_fastq = $datapath."Processed/".$key."/".$`.".fastq";
            }
         else{die "Something is wrong!017\n";}
      if($R1files{$key} =~ m/_R1_/s)
         {
           $outputfile2 = "$datapath"."Processed/".$key."/".$`;
           $log_file2 = "$datapath"."Processed/".$key."/".$`."_log";
         }
        else{die "Something is wrong!018\n";}
    # Forks and returns the pid for the child:
         my $pid = $pm->start and next;
    # we are now in the child process
         system $code, "--I1", $I1_fastq, "--R2", $R2_fastq, "--out", $outputfile2, "--bar", $bar_file, "--end-I1", $I1_end, "--end-R2", $R2_end, "--primer-I1", $I1_length, "--primer-R2", $R2_length, "--middle", $keep_middle, "--pip", $pippath, "--soft", $nucleotide, "--log", $log_file2;
         $pm->finish; # Terminates the child process
       }
   }
$pm->wait_all_children;
##Finish BarBIQ_M_I1R2.pl##

print "Barcode_Step_1.pl has done at:";
print scalar localtime;
print "\n";


##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.2.0
#2022.10.13
