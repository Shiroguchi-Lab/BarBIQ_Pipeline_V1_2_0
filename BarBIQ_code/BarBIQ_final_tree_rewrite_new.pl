#! /usr/bin/env perl
### sort the clusterizer file by the illumina read name

use strict;
use warnings;

my $log_file = $ARGV[1];
open (my $LOG, '>>', $log_file);
select $LOG;

print "Now you are runing $0\n";
print "The parameters are: @ARGV\n";

my $inputfile = "$ARGV[0]"; ## The sum clusterizer result file
my $outputfile = "$ARGV[0]_Rewrite"; ##Sorted file 
unlink $outputfile;

open(FILE, $inputfile) or die "cannot open input file '$inputfile' $!";
open(OUTF, '>>', $outputfile) or die "cannot open input file '$outputfile' $!";
print OUTF ("ID\tdomain\tdomain_s\tphylum\tphylum_s\tclass\tclass_s\torder\torder_s\tfamily\tfamily_s\tgenus\tgenus_s\tspecies\tspecies_s\n");
while(my $gi=<FILE>)
     {
      chomp $gi;
      my @info=split(/\s+/,$gi);
      for(my $i=0; $i<=$#info; $i++)
        {
          if($info[$i] =~ /\A"/){
                    if($' =~ /"\Z/){
                         $info[$i] = $`;
                      }
                  }
        }
      my @taxa_IP =  &taxa(\@info);
      &rewrite(\@info, \@taxa_IP);
     }

sub taxa {
     my ($data) = $_[0];
     my $domain = 0;
     my $phylum = 0;
     my $class = 0;
     my $order = 0;
     my $family= 0;
     my $genus = 0;
     my $subclass =0; 
     my $suborder =0;
     for (my $i=0; $i<=$#{$data}; $i++)
        {
          if($$data[$i] eq "domain") { $domain = $i;} 
         elsif ($$data[$i] eq "phylum") { $phylum = $i;}
         elsif ($$data[$i] eq "class") { $class = $i;}
         elsif ($$data[$i] eq "order") { $order = $i;}
         elsif ($$data[$i] eq "family") { $family = $i;}
         elsif ($$data[$i] eq "genus") { $genus = $i;}
         elsif ($$data[$i] eq "subclass") { $subclass = $i;}
         elsif ($$data[$i] eq "suborder") { $suborder = $i;}
        }
     return ($domain,$phylum,$class,$order,$family,$genus, $subclass, $suborder);
    }

sub rewrite {
     my ($data, $ID) = @_;
     print OUTF ($$data[0]);
     my $start = 5;
     # domian level
     if ($$ID[0]) { print OUTF ("\t$$data[$$ID[0]-1]\t$$data[$$ID[0]+1]"); 
                  $start = $$ID[0]+2; }
    else { print "$$data[0] does not contain the domain name\n"; 
           print OUTF ("\t*\t-"); }
     #phylum level
     if ($$ID[1]) { 
            print OUTF ("\t");
            my $head=0;
                for( my $i= $start; $i<$$ID[1]; $i++)
                         {
                           if($head) { print OUTF ("_");$head = 1; }
                           print OUTF ("$$data[$i]");
                         }
            print OUTF ("\t$$data[$$ID[1]+1]");       
                  $start = $$ID[1]+2; }
     else {
           print "$$data[0] does not contain the phylum name\n";
           print OUTF ("\t*\t-");
          }
     #class level
     if ($$ID[2]) { 
           print OUTF ("\t");
           my $head=0;
               for( my $i= $start; $i<$$ID[2]; $i++)
                         { 
                           if($head) { print OUTF ("_");$head = 1; }
                           print OUTF ("$$data[$i]");
                         }
                  print OUTF ("\t$$data[$$ID[2]+1]");
                  $start = $$ID[2]+2; }
     else {
           print "$$data[0] does not contain the class name\n";
           print OUTF ("\t*\t-");
          }
     #order level
     if($$ID[6]) {print "$$data[0] contain subclass name\n";
                  $start = $$ID[6] + 2;
                 }
     if ($$ID[3]) {
           print OUTF ("\t");
           my $head=0;
               for( my $i= $start; $i<$$ID[3]; $i++)
                         { 
                          if($head) { print OUTF ("_");$head = 1; } print OUTF ("$$data[$i]");
                         }
                  print OUTF ("\t$$data[$$ID[3]+1]");
                  $start = $$ID[3]+2; }
     else {
           print "$$data[0] does not contain the order name\n";
           print OUTF ("\t*\t-");
          }
     #family level
     if($$ID[7]) {print "$$data[0] contain suborder name\n";
                  $start = $$ID[7] + 2;
                 }
     if ($$ID[4]) {
           print OUTF ("\t");
           my $head=0;
               for( my $i= $start; $i<$$ID[4]; $i++)
                         {
                           if($head) { print OUTF ("_");$head = 1; }
                           print OUTF ("$$data[$i]");
                         }
                  print OUTF ("\t$$data[$$ID[4]+1]");
                  $start = $$ID[4]+2; }
     else {
           print "$$data[0] does not contain the family name\n";
           print OUTF ("\t*\t-");
          }
     #genus level
     if ($$ID[5]) {
           print OUTF ("\t");
           my $head=0;
               for( my $i= $start; $i<$$ID[5]; $i++)
                         {
                           if($head) { print OUTF ("_");$head = 1; }
                           print OUTF ("$$data[$i]");
                         }
                  print OUTF ("\t$$data[$$ID[5]+1]");
                  $start = $$ID[5]+2; }
     else {
           print "$$data[0] does not contain the genus name\n";
           print OUTF ("\t*\t-");
          }
     print OUTF ("\t$$data[0]\t0\n");
    }

close OUTF;
close FILE;

##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.20
