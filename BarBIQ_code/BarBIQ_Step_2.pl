#! /usr/bin/env perl
##########################################################################################################
#####Description of this code#####
# This code is the second step of the BarBIQ pipline.
# It is a combined pipeline of  
#    BarBIQ_final_repeat_no_junk_deletion.pl
#    BarBIQ_final_repeat_identity_Editdis.pl
#    BarBIQ_final_repeat_delete_low_count.pl 
#    BarBIQ_final_merge_all_samples.pl    
#    BarBIQ_final_library.pl 
#    BarBIQ_final_retrieving_false_negative_sequence_types.pl
#    BarBIQ_final_droplets_overlap.pl
#    BarBIQ_final_merge_all_overlaps.pl
#    BarBIQ_final_fitting_OD.r
#    BarBIQ_final_overlap_selete_by_PV_seq.pl    
#    BarBIQ_final_overlap_selete_by_PV_seq_threshold.pl
#    BarBIQ_final_overlap_selete_by_PV_step3.pl
#    BarBIQ_final_overlap_groups.pl
#    BarBIQ_final_lib_COTU_ID.pl
#    BarBIQ_final_groups_COTU_ID.pl
#    BarBIQ_final_add_bacteria_name_to_final_rep_seq_file.pl
#    BarBIQ_final_count_bacteria.pl
#    BarBIQ_final_compare_bacteria_count.pl
#    BarBIQ_final_review_contamination_bac_123rep.pl
#    BarBIQ_final_compare_datasets.pl
#    BarBIQ_final_lib_COTU_clean.pl
#    BarBIQ_final_LIB_fastafile.pl
#    BarBIQ_final_tree_rewrite_new.pl
#    BarBIQ_final_Taxonomy_groups.pl
#    BarBIQ_final_Taxonomy_COTU.pl 
#    BarBIQ_final_add_Taxonomy_to_bac_count_publish.pl
#    BarBIQ_final_bac_con_nomalization_xlsx.pl
#    BarBIQ_final_bar_seq_xlsx.pl 
# So please make sure these codes are all under the same directory of this combined code.
# If you use this combined code instead of those separated codes respectively, meaning you want to use all default parameters without any modification.
# The default parameters please check those codes in detail respectively. 
# This code uses these input files(see examples), please prepare in advance:
# BarBIQ_example_inputfile_2.txt: information about all samples you are going to compare
# BarBIQ_example_parameter_2.txt : necessary parameters
###########################################################################################################
#####how to run this code #####
##command##
## BarBIQ_Step_2.pl --in BarBIQ_example_inputfile_2.txt --p BarBIQ_example_parameter_2.txt (--middle Yes)
##interpretation##
# --in: a file should be prepared like the example BarBIQ_example_inputfile_2.txt 
# --p: a file including necessary parameters should be prepared like the example BarBIQ_example_parameter_2.txt
# --middle: if you want to keep the middlefiles which generated during processing, please set this parameter as "Yes", default is "No"
########################################################################################################## 

use strict;
use warnings;
use IPC::System::Simple qw(system);
use File::Path qw(make_path);
use File::Copy;
use Parallel::ForkManager;


##read command##
my ($i,$inputfile,$parefile);
my $keep_middle="No";
for($i=0; $i<=$#ARGV; $i=$i+2)
    {
        if ($ARGV[$i] eq "--in")  {$inputfile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--p") {$parefile = $ARGV[$i+1];}
     elsif ($ARGV[$i] eq "--middle") {$keep_middle = $ARGV[$i+1];}
     else                         {die "Your input is wrong!!!\n Please input \"--in BarBIQ_example_inputfile_2.txt --p BarBIQ_example_parameter_2.txt --out outputfile\"\n $!";}
    }
if(!$inputfile)   {die "Your input is wrong!!!\n Please input \"--in: BarBIQ_example_inputfile_2.txt\"\n $!";}
if(!$parefile)  {die "Your input is wrong!!!\n Please input \"--p BarBIQ_example_parameter_2.txt\"\n $!";}
##read command##

##Check the parameter file##
my ($MAX_PROCESSES,$Output_path);
my ($gi, @info, $pippath, $code, $RDPpath, $RDP_c);
open(FILE,$parefile) or die "Could not open file '$parefile' $!"; # open inputfile
while(my $gi=<FILE>)  ## read Real-data file names
    {
     chomp $gi;
     my @info=split(/\s+/,$gi);
     if($info[0] eq "Perl_code_path:") {$pippath = $info[1];}
  elsif($info[0] eq "Output_path:") {$Output_path = $info[1];}
  elsif($info[0] eq "RDPTools:") {$RDPpath = $info[1];}
  elsif($info[0] eq "RDPTools_c:") {$RDP_c = $info[1];}
  elsif($info[0] eq "Max_processes:") {$MAX_PROCESSES = $info[1];}
  elsif($info[0] ne "#")
      {
        print ("$info[0]\n");
      }
  }
close FILE;
if(!$Output_path) {die "Your parameter file is wrong!!\n Please privide Output_path\n $!"};
if (-e $Output_path && -d $Output_path) {
    print "The output folder $Output_path is existed\n";
   }
else{
    make_path($Output_path);
    print "The output folder $Output_path is newly made\n";
    }
my @parts = split('/', $Output_path);
my $outputfile1=$parts[$#parts];
my $log_file = $Output_path.$outputfile1."_log";
unlink $log_file;
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are running program: $0\n";
print "Your parameters are: @ARGV\n";
print "Started at: ";
print scalar localtime;
print "\n";

if(!$pippath)  {die "Your parameter file is wrong!!\n Please privide Perl_code_path\n $!";}
if(!$MAX_PROCESSES) {die "Your parameter file is wrong!!\n Please privide MAX_PROCESSES\n $!";}
##Finish check the parameter file##


##check the inputfile##
my (%folder, %index, %SampleID, %Repeat, %Sampletypes,%MutiSeqDetect,%Control);
# my (@real_data, %index, @std_data, $fixed_barcode_file, %samplename, %R2files, %I1files, %R1files, %Sampletypes, $fixed_barcode_std_file, $gi, @info, $pippath, $datapath, $code);
open(FILE,$inputfile) or die "Could not open file '$inputfile' $!"; # open inputfile
chomp ($gi=<FILE>);
@info=split(/\s+/,$gi);
my ($foldcol, $indexcol, $Typecol,$Samplecol, $Repeatcol,$Muticol,$Concol,$TTcol,$SDcol);
for($i=0; $i<=$#info; $i++)
     {  
        if($info[$i] eq "Index"){$indexcol = $i;}
     elsif($info[$i] eq "Folder"){$foldcol = $i;}
     elsif($info[$i] eq "SampleID"){$Samplecol = $i;}
     elsif($info[$i] eq "Repeat"){$Repeatcol = $i;}
     elsif($info[$i] eq "SampleType"){$Typecol = $i;}
     elsif($info[$i] eq "Multi-seq-detect"){$Muticol = $i;}
     elsif($info[$i] eq "Control"){$Concol = $i;}
     elsif($info[$i] eq "Total-conc(copies/mg)"){$TTcol = $i;}
     elsif($info[$i] eq "Total-conc-SD"){$SDcol = $i;}
     #else{die "Your input is wrong!!!\n The BarBIQ_example_inputfile.txt should have these columns:Index	I1_file	R1_file	R2_file	SampleType \n $!";}
    }
if(defined $foldcol && defined $indexcol && defined $Typecol && defined $Samplecol && defined $Repeatcol && $Muticol && $Concol && $TTcol && $SDcol) 
      {
      }
else {
         die "Your input is wrong!!!\n The BarBIQ_example_inputfile.txt should have these columns:Folder  Index   SampleType      SampleID        Repeat	Multi-seq-detect        Control\n $!";
      }
my $order = 0;
my $total_conc = $Output_path."Total_concentration.txt";
open(OUTF, '>>', $total_conc) or die "Could not open file '$total_conc' $!";
print OUTF ("Sample\tConc.(copies/mg)\tS.d\n");
while($gi=<FILE>)  ## read Real-data file names
    {
     $order++;
     chomp $gi;
     @info=split(/\s+/,$gi); 
     if($info[$Typecol] ne "SpikeIn")
        { 
       $folder{$order} = $info[$foldcol];
       $index{$order}=$info[$indexcol];
       $SampleID{$order}=$info[$Samplecol];
       $Repeat{$order}=$info[$Repeatcol];
       $Sampletypes{$order}=$info[$Typecol];
       $MutiSeqDetect{$order}=$info[$Muticol];
       $Control{$order}=$info[$Concol];
       if($info[$Typecol] ne "Blank")
          { print OUTF ("$info[$Samplecol]\t$info[$TTcol]\t$info[$SDcol]\n"); }
        }
    }
close FILE; # close inputfile
close OUTF;
if(!%index) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%folder) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%SampleID) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%Repeat) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
if(!%Sampletypes) {die "Your inputfile '$inputfile' is wrong!!! please check!!! $!";}
print "Inputfile is OK!\n";
##Finish check the inputfile##


##Run BarBIQ_final_repeat_no_junk_deletion.pl##
print "###############################################################################\n";
my @unique;
my %seen;
my %outputname2;
$code="$pippath"."BarBIQ_final_repeat_no_junk_deletion.pl"; 
foreach my $key (sort keys %SampleID) 
     {
       if (! $seen{$SampleID{$key}}++ && ($Sampletypes{$key} ne "Blank")) {
       push @unique, $SampleID{$key};
      }
     }
#print @unique;
for($i=0; $i<=$#unique; $i++)
     {
        print "Now you are processing the sample $unique[$i]\n";
        $outputname2{$unique[$i]} = $Output_path.$unique[$i]."_step12";
        my @ppp = $outputname2{$unique[$i]};
        print "Your output file is $outputname2{$unique[$i]}\n";
        foreach my $key2 (sort keys %SampleID)
            {
                if($SampleID{$key2} eq $unique[$i])
                   {
                     my $datapath = $folder{$key2}."Processed/".$index{$key2}."/";
                     opendir my $dir, $datapath or die "Cannot open directory: $datapath $!";
                     my @allfiles = readdir $dir;
                     closedir $dir;
                     my @array = grep {/_clean_sub1_sub2_shift_link_derr_indels_chimeras_similar_statistic\z/} @allfiles;
                     my $file = $folder{$key2}."Processed/".$index{$key2}."/".$array[0];
                     print "Your inputfile is $file\n";
                     unshift(@ppp, $file);
                     #my $ppp2 = join(',', @ppp);
                     #print "$ppp2\n";
                    }
             }
         system $code,@ppp,$log_file;
      }

##Run BarBIQ_final_repeat_identity_Editdis.pl and BarBIQ_final_repeat_delete_low_count.pl##
$code="$pippath"."BarBIQ_final_repeat_identity_Editdis.pl";
my $code2 = "$pippath"."BarBIQ_final_repeat_delete_low_count.pl"; 
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);
foreach my $key (keys %outputname2)
     {
        my $pid = $pm->start and next;
        my $inputfile=$outputname2{$key}."_normalization";
        system $code, $inputfile,$log_file;
        $inputfile=$inputfile."_clean";
        system $code2, $inputfile, $log_file;
        $pm->finish; # Terminates the child process
     }
$pm->wait_all_children;

if($keep_middle eq "No")
     {
       foreach my $key (keys %outputname2)
         {
            unlink $outputname2{$key}."_normalization";
            unlink $outputname2{$key}."_normalization_clean";
            unlink $outputname2{$key}."_normalization_identity";
         }
     }
##Finished BarBIQ_final_repeat_identity_Editdis.pl and BarBIQ_final_repeat_delete_low_count.pl##

##Run BarBIQ_final_merge_all_samples.pl##
my @samples;
foreach my $key (keys %outputname2)
     {
       my $file=$outputname2{$key}."_normalization_clean_del_low";
       unshift(@samples, $file);
     }
my $outputfile3 = $Output_path."All_samples_step13";
$code = "$pippath"."BarBIQ_final_merge_all_samples.pl";
system $code, @samples, $outputfile3, $log_file;
if($keep_middle eq "No")
     {
       foreach my $key (keys %outputname2)
         {
            unlink $outputname2{$key}."_normalization_clean_del_low";
         }
     }
##Finished BarBIQ_final_merge_all_samples.pl## 

##Run BarBIQ_final_library.pl##
$code = "$pippath"."BarBIQ_final_library.pl";
my $outputfile4 = $Output_path."All_samples_BarSeqs";
system $code, $outputfile3, $outputfile4,$log_file;
if($keep_middle eq "No")
     {
       unlink $outputfile3;
     }
##Finished BarBIQ_final_library.pl##

##Run BarBIQ_final_retrieving_false_negative_sequence_types.pl##
my @unique2;
foreach my $key (sort keys %SampleID)
     {
       if ($MutiSeqDetect{$key} eq "Yes" || $MutiSeqDetect{$key} eq "yes") {
       push @unique2, $key;
       }
     }
$code = "$pippath"."BarBIQ_final_retrieving_false_negative_sequence_types.pl";
$code2 = "$pippath"."BarBIQ_final_droplets_overlap.pl";
$pm = Parallel::ForkManager->new($MAX_PROCESSES);
for($i=0; $i<=$#unique2; $i++)
     {
        my $pid = $pm->start and next;
        my $datapath = $folder{$unique2[$i]}."Processed/".$index{$unique2[$i]}."/";
        opendir my $dir, $datapath or die "Cannot open directory: $datapath $!";
        my @allfiles = readdir $dir;
        closedir $dir;
        my @array = grep {/_clean_sub1_0.1_sub2_shift_link\z/} @allfiles;
        my $file = $folder{$unique2[$i]}."Processed/".$index{$unique2[$i]}."/".$array[0];
        print "Your inputfile is $file\n";
        system $code, '--in', $file, '--bar', $outputfile4,'--log',$log_file;
        my $file2 = $file."_Rev";
        system $code2, '--in', $file2, '--bar', $outputfile4,'--log',$log_file;
        if($keep_middle eq "No")
           {
             unlink $file2;
           }
        $pm->finish; # Terminates the child process
     }
$pm->wait_all_children;

##Finished BarBIQ_final_retrieving_false_negative_sequence_types.pl## 

##Run BarBIQ_final_merge_all_overlaps.pl##
print "###############################################################################\n";
my @samples2;
for($i=0; $i<=$#unique2; $i++)
     {
        my $datapath = $folder{$unique2[$i]}."Processed/".$index{$unique2[$i]}."/";
        opendir my $dir, $datapath or die "Cannot open directory: $datapath $!";
        my @allfiles = readdir $dir;
        closedir $dir;
        my @array = grep {/_clean_sub1_0.1_sub2_shift_link_Rev_DropOver\z/} @allfiles;
        my $file = $folder{$unique2[$i]}."Processed/".$index{$unique2[$i]}."/".$array[0];
        my $file2 = $Output_path."Sample_S".$unique2[$i]."_Rev_DropOver";
        copy($file, $file2);        
        push(@samples2, $file2);
     }
my $outputfile5 = $Output_path."All_samples_Overlap";
$code = "$pippath"."BarBIQ_final_merge_all_overlaps.pl";
system $code, @samples2, $outputfile5, $log_file;
if($keep_middle eq "No")
           {
            foreach my $file (@samples2) {
             unlink $file;
             }
           }
##Finished BarBIQ_final_merge_all_overlaps.pl##

##Run BarBIQ_final_fitting_OD.r##
my $outputfile6 = $Output_path."Index_Overlap";
unlink $outputfile6;
open(OUTF, '>>', $outputfile6) or die "cannot open input file '$outputfile6' $!";
my @indexes;
for($i=0; $i<=$#unique; $i++)
     {
       my $ins = "S".$unique2[$i];
       push @indexes, $ins;
     }
my $in = join("\t", @indexes);
print OUTF ("$in\n");
close OUTF;
my $outputfile7 = $Output_path."EDrops.txt";
$code = "$pippath"."BarBIQ_final_fitting_OD.r";
system "Rscript", "--vanilla", $code, $outputfile5, $outputfile6, $outputfile7,$log_file;
if($keep_middle eq "No")
           {
             unlink $outputfile6;
           }
##Finished BarBIQ_final_fitting_OD.r##

##Run BarBIQ_final_overlap_selete_by_PV_seq.pl##
$code = "$pippath"."BarBIQ_final_overlap_selete_by_PV_seq.pl";
my $pv = "$pippath"."Estimated_confidence_intervals_of_log10_Poisson_Overlap_by_simulation";
system $code, '--overlap', $outputfile5, '--PV', $pv, '--EDrop', $outputfile7, '--log', $log_file;
##Finished BarBIQ_final_overlap_selete_by_PV_seq.pl##

##Run BarBIQ_final_overlap_selete_by_PV_seq_threshold.pl##
my $outputfile8="$outputfile5"."_PV_list";
$code = "$pippath"."BarBIQ_final_overlap_selete_by_PV_seq_threshold.pl";
my $threshold = 0.5; # How much portion of samples show the Overlap is higher than Poission distribution
system $code, $outputfile8, $threshold, $log_file;
##Finished BarBIQ_final_overlap_selete_by_PV_seq_threshold.pl##

##Run BarBIQ_final_overlap_selete_by_PV_step3.pl##
my $outputfile9="$outputfile8"."_$threshold";
$code = "$pippath"."BarBIQ_final_overlap_selete_by_PV_step3.pl";
system $code, '--overlap', $outputfile5, '--select', $outputfile9, '--log', $log_file;
##Finished BarBIQ_final_overlap_selete_by_PV_step3.pl##

##Run BarBIQ_final_overlap_groups.pl##
my $outputfile10="$outputfile5"."_select";
$code = "$pippath"."BarBIQ_final_overlap_groups.pl";
system $code, $outputfile10, $log_file;
##Finished BarBIQ_final_overlap_groups.pl##

##Run BarBIQ_final_lib_COTU_ID.pl##
my $outputfile11="$outputfile10"."_groups";
$code = "$pippath"."BarBIQ_final_lib_COTU_ID.pl";
system $code, '--bar', $outputfile4, '--group', $outputfile11, '--log', $log_file;
##Finished BarBIQ_final_lib_COTU_ID.pl##

##Run BarBIQ_final_groups_COTU_ID.pl##
my $outputfile12="$outputfile4"."_COTU";
$code = "$pippath"."BarBIQ_final_groups_COTU_ID.pl";
system $code, '--bar', $outputfile12, '--group', $outputfile11, '--log', $log_file;
my $outputfile112=$outputfile11."_COTU";
##BarBIQ_final_groups_COTU_ID.pl##

if($keep_middle eq "No")
           {
             unlink ($outputfile4,$outputfile5,$outputfile8,$outputfile9,$outputfile10,$outputfile11);
           }


##Run BarBIQ_final_retrieving_false_negative_sequence_types.pl##
##Run BarBIQ_final_add_bacteria_name_to_final_rep_seq_file.pl##
##Run BarBIQ_final_count_bacteria.pl##
print "###############################################################################\n";
my @unique3 = sort keys %SampleID;
# print @unique3;
$code = "$pippath"."BarBIQ_final_retrieving_false_negative_sequence_types.pl";
$code2 = "$pippath"."BarBIQ_final_add_bacteria_name_to_final_rep_seq_file.pl";
my $code3 = "$pippath"."BarBIQ_final_count_bacteria.pl";
$pm = Parallel::ForkManager->new($MAX_PROCESSES);
for($i=0; $i<=$#unique3; $i++)
     {
        my $pid = $pm->start and next;
        my $datapath = $folder{$unique3[$i]}."Processed/".$index{$unique3[$i]}."/";
        opendir my $dir, $datapath or die "Cannot open directory: $datapath $!";
        my @allfiles = readdir $dir;
        closedir $dir;
        my @array = grep {/_clean_sub1_sub2_shift_link_derr_indels_chimeras_similar\z/} @allfiles;
        my $file = $folder{$unique3[$i]}."Processed/".$index{$unique3[$i]}."/".$array[0];
        print "Your inputfile is $file\n";
        system $code, '--in', $file, '--bar', $outputfile12,'--log',$log_file;
        my $file2 = $file."_Rev";
        system $code2, '--in', $file2, '--COTU', $outputfile12,'--log',$log_file;
        my $file3="$file2"."_ID";
        system $code3, '--in', $file3,'--log',$log_file;  
        if($keep_middle eq "No")
           {
             unlink $file2;
             unlink $file3;
           }
        $pm->finish; # Terminates the child process
     }
$pm->wait_all_children;
##Finished BarBIQ_final_add_bacteria_name_to_final_rep_seq_file.pl##
##Finished BarBIQ_final_count_bacteria.pl##
##Finished BarBIQ_final_retrieving_false_negative_sequence_types.pl## 

##Run BarBIQ_final_compare_bacteria_count.pl##
##Run BarBIQ_final_review_contamination_bac_123rep.pl##
print "###############################################################################\n";
my @unique4;
foreach my $key (sort keys %Control)
     {
       if (! $seen{$Control{$key}}++) {
             push @unique4, $Control{$key};
           }
     }
$code = "$pippath"."BarBIQ_final_compare_bacteria_count.pl";
$code2 = "$pippath"."BarBIQ_final_review_contamination_bac_123rep.pl";
my %outputname17;
# print @unique4;
for($i=0; $i<=$#unique4; $i++)
    {
      print "Your are processing the control group $unique4[$i]\n";
      $outputname17{$unique4[$i]} = $Output_path."Cell_count_".$unique4[$i];
      my @ppp = $outputname17{$unique4[$i]};
      print "Your output file is $outputname17{$unique4[$i]}\n";
      my @file2;
      my $outputfile172 = $outputname17{$unique4[$i]}."_IF";
      unlink $outputfile172;
      open(OUTF, '>>', $outputfile172) or die "Could not open file '$outputfile172' $!"; 
      print OUTF ("Index\tSample\tType\n");
      foreach my $key2 (reverse sort keys %Control)
            {
                if($Control{$key2} eq $unique4[$i])
                   {
                     my $datapath = $folder{$key2}."Processed/".$index{$key2}."/";
                     opendir my $dir, $datapath or die "Cannot open directory: $datapath $!";
                     my @allfiles = readdir $dir;
                     closedir $dir;
                     my @array = grep {/_clean_sub1_sub2_shift_link_derr_indels_chimeras_similar_Rev_ID_bac_count\z/} @allfiles;
                     my $file = $folder{$key2}."Processed/".$index{$key2}."/".$array[0];
                     my $file2 = $Output_path."Sample_S".$key2."_Bac_Count";
                     copy($file, $file2);
                     print "Your inputfile is $file\n";
                     if($Sampletypes{$key2} eq "Blank")
                         {
                            print OUTF ("S$key2\t$SampleID{$key2}\t$Sampletypes{$key2}\n");
                         }
                     else{print OUTF ("S$key2\t$SampleID{$key2}\t$SampleID{$key2}\n");}
                     unshift(@ppp, $file2);
                     push(@file2, $file2);
                     #my $ppp2 = join(',', @ppp);
                     #print "$ppp2\n";
                    }
             }
     close OUTF;
     system $code,@ppp,$log_file; 
     if($keep_middle eq "No")
           {
            foreach my $file (@file2) {
             unlink $file;
             }
           }
     system $code2,'--in',$outputname17{$unique4[$i]}, '--IF', $outputfile172, '--log', $log_file;     
      }
##Finished BarBIQ_final_compare_bacteria_count.pl##
##Finished BarBIQ_final_review_contamination_bac_123rep.pl##

##Run BarBIQ_final_compare_datasets.pl##
$code = "$pippath"."BarBIQ_final_compare_datasets.pl";
my $outputfile173 = $Output_path."All_samples_cell_count";
my @ppp = $outputfile173;
foreach my $key3 (sort keys %outputname17)
      {
        unshift(@ppp, $outputname17{$key3}."_ave_del_contamination");
      }
system $code, @ppp, $log_file;
##Finished BarBIQ_final_compare_datasets.pl##

##Run BarBIQ_final_lib_COTU_clean.pl##
$code = "$pippath"."BarBIQ_final_lib_COTU_clean.pl";
system $code,'--bcc', $outputfile173, '--COTULIB', $outputfile12, '--log', $log_file;
##Finished BarBIQ_final_lib_COTU_clean.pl##

##Run BarBIQ_final_LIB_fastafile.pl##
$code = "$pippath"."BarBIQ_final_LIB_fastafile.pl";
my $outputfile13 = $outputfile12."_clean";
system $code,'--lib', $outputfile13, '--log', $log_file;
##Finished BarBIQ_final_LIB_fastafile.pl##

if($keep_middle eq "No")
     {
       foreach my $file (keys %outputname17) {
       unlink $outputname17{$file};
       unlink $outputname17{$file}."_ave_del_contamination";
       unlink $outputname17{$file}."_IF";
     }
   }


##Run RDP classifier##
my $outputfile14 = $outputfile13.".fa";
my $outputfile15 = $outputfile13."_taxanomy.txt";
print "You are using assignment confidence cutoff $RDP_c for the RDP classifier\n";
$code = "$RDPpath"."classifier.jar";
system 'java', '-jar', $code, 'classify', '-c', $RDP_c, '-o', $outputfile15, $outputfile14; 
$code = "$pippath"."BarBIQ_final_tree_rewrite_new.pl";
system $code,$outputfile15,$log_file;
my $outputfile16 = $outputfile15."_Rewrite";
$code = "$pippath"."BarBIQ_final_Taxonomy_groups.pl";
system $code,'--group', $outputfile112,'--taxa',$outputfile16, '--log', $log_file;
my $outputfile17 = $outputfile112."_RDPpdt_taxonomy";
$code = "$pippath"."BarBIQ_final_Taxonomy_COTU.pl";
system $code,'--lib', $outputfile12,'--group', $outputfile17,'--taxa',$outputfile16, '--log', $log_file;
##Finished RDP classifier##

##Run BarBIQ_final_add_Taxonomy_to_bac_count_publish.pl##
$code = "$pippath"."BarBIQ_final_add_Taxonomy_to_bac_count_publish.pl";
my $outputfile18 = $outputfile16."_COTU";
system $code,'--bacc', $outputfile173,'--taxa',$outputfile18, '--log', $log_file;
##Finished BarBIQ_final_add_Taxonomy_to_bac_count_publish.pl##

if($keep_middle eq "No")
     {
       unlink ($outputfile14, $outputfile15, $outputfile16, $outputfile18);
     }

##Run BarBIQ_final_bac_con_nomalization_xlsx.pl##
my $outputfile20 = $outputfile173."_annotation_RDP_classifier.txt";
my $outputfile19 = $Output_path."BarBIQ_Quantification.xlsx";
$code = "$pippath"."BarBIQ_final_bac_con_nomalization_xlsx.pl";
system $code, '--count', $outputfile20, '--total', $total_conc, '--out', $outputfile19, '--taxa', 'Classifier', '--cutoff', $RDP_c, '--log', $log_file;
##Finished BarBIQ_final_bac_con_nomalization_xlsx.pl ##

##Run BarBIQ_final_bar_seq_xlsx.pl##
my $outputfile21 = $Output_path."BarBIQ_Bar_sequences.xlsx";
$code = "$pippath"."BarBIQ_final_bar_seq_xlsx.pl";
system $code, '--bar', $outputfile13, '--out', $outputfile21,'--log', $log_file;
##Finished BarBIQ_final_bar_seq_xlsx.pl##

##Run BarBIQ_final_bar_seq_fasta_with_taxanomy.pl##
my $outputfile22 = $Output_path."BarBIQ_Bar_sequences.fasta";
$code = "$pippath"."BarBIQ_final_bar_seq_fasta_with_taxanomy.pl";
system $code, '--bar', $outputfile13, '--taxa', $outputfile20, '--out', $outputfile22,'--log', $log_file;
##Finished BarBIQ_final_bar_seq_fasta_with_taxanomy.pl##

print "BarBIQ_Step_2.pl has done at:";
print scalar localtime;
print "\n";

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.2.2
#2023.02.15
