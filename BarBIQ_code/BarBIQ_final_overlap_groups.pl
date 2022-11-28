#! /usr/bin/env perl
#####Description of this code#####
#This code is used to find mutipule Bar sequences from the same bacterium, if we found common Bar sequences in different same-bacterial Bar sequences, we concluded that all these Bar sequences were in one bacterium. 
################################################################################################################how to run this code #####
##command###
#BarBIQ_final_overlap_groups.pl inputfile
##explaination##
#inputfile: output file from BarBIQ_final_overlap_selete_by_PV_step3.pl
########################################################################################################## 
#####Install#####
##None
############################################################################################################
#####code##### 

use strict;
use warnings;

my $log_file = $ARGV[1];
open (my $LOG, '>>', $log_file);
select $LOG;

print "now you are running program: $0\n";
print "Your parameters are: @ARGV\n";

my $inputfile = "$ARGV[0]"; ## 
my $outputfile = "$ARGV[0]"."_groups"; ##Sorted file 
my $outputfile_nodes = "$ARGV[0]"."_nodes.csv";
my $outputfile_edges = "$ARGV[0]"."_edges.csv";
unlink $outputfile;
unlink $outputfile_nodes;
unlink $outputfile_edges;
my %group;
open(IN, $inputfile) || die "canot open input file '$inputfile' $!";
open(OUTF, '>>', $outputfile) or die "Could not open file '$outputfile' $!";
my $gi = <IN>;
chomp $gi;
my @info=split(/\s+/,$gi);
my $ratio;
my $SEQ1;
my $SEQ2;
for (my $i=0; $i<=$#info; $i++)
   {
     if($info[$i] eq "SEQ1") {$SEQ1 = $i;}
     if($info[$i] eq "Ratio") {$ratio = $i;}
     if($info[$i] eq "SEQ2") {$SEQ2 = $i;}
   }
if (!($SEQ1 && $ratio && $SEQ2)) {die "File '$inputfile' is wrong!!!001\n";}
my $i=1;
my %nodes;
my @edges;
push @edges, [("from","to","ratio")];
while($gi = <IN>)
   {
    #   print "$i\n";
     chomp $gi;
     @info=split(/\s+/,$gi);
     push @edges, [($info[0],$info[2],$info[$ratio])];
     $nodes{$info[0]} = $info[$SEQ1];
     $nodes{$info[2]} = $info[$SEQ2];
     if(%group)
       {
         my $add=1;
         foreach my $key (keys %group)
              {
                if($add)
                 {
                  my @a=split(/\s+/,$group{$key});
                  my @b=($info[0], $info[2]);
                  my %hash_a = map{$_=>1} @a;
                  my %hash_b = map{$_=>1} @b;
                  my %merge_all = map {$_ => 1} @a,@b;
                  my @common = grep {$hash_a{$_}} @b;
                  my @merge = keys (%merge_all);
                  if($common[0]) { my $m=join("\t", @merge); $group{$key} = $m;$add=0;}
                 }
              }
          if($add) {$group{$i}="$info[0]\t$info[2]"; $i++;}
       }     
     else{$group{$i}="$info[0]\t$info[2]"; $i++;}
#print "$i\n";
   }

my $overlap=1;
use Data::Dumper;
while($overlap)
{
$overlap=0;
foreach my $key1 (keys %group)
   {
      foreach my $key2 (keys %group)
         {
           if ($key1 lt $key2)
               {
                 my @a=split(/\s+/,$group{$key1});
                 my @b=split(/\s+/,$group{$key2});
                 my %hash_a = map{$_=>1} @a;
                 my %hash_b = map{$_=>1} @b;
                 my %merge_all = map {$_ => 1} @a,@b;
                 my @common = grep {$hash_a{$_}} @b;
                 my @merge = keys (%merge_all);
                 if($common[0]) 
                     {
                      my $m=join("\t", @merge);
                      $group{$key1} = $m;
                      delete $group{$key2};
                      $overlap=1;
                      print "Common :\n";
                      print Dumper(\@common);
                      print "A :\n";
                      print Dumper(\@a);
                      print  "B :\n";
                      print Dumper(\@b);
                      last;
                     }
               }    
         }
     if($overlap) {last;} 
   } 
} 
my $p=1;
foreach my $key (sort keys %group)
    {
     print OUTF ("Group_$p\t$group{$key}\n");
     $p++;
    }
close OUTF;
close IN;

open(OUTFE, '>>', $outputfile_edges) or die "Could not open file '$outputfile_edges' $!";
for my $row ( @edges ) {
    print OUTFE "\"";
    print OUTFE join( '","', @$row );
    print OUTFE "\"\n";
}
close OUTFE;

open(OUTFN, '>>', $outputfile_nodes) or die "Could not open file '$outputfile_nodes' $!";
print OUTFN ("\"LIBID\",\"Sequence\"\n");
foreach my $key (sort keys %nodes)
   {
     print OUTFN ("\"$key\",\"$nodes{$key}\"\n");
   }
close OUTFN;

##end##

#####Author#####
#Jianshi Frank Jin

#####Version#####
#V1.001
#2018.07.27
