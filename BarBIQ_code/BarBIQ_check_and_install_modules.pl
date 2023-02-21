#! /usr/bin/env perl
#################################################################################################
#####Description of this code#####
#This code is used to install all the required modules for BarBIQ pipeline.
#################################################################################################
#####how to run this code#####
##command##
#BarBIQ_check_and_install_modules.pl
####################################################################################################

#####code#######
use strict;
use warnings;
use CPAN;

my @modules = ("Bundle::BioPerl", "IPC::System::Simple", "File::Path", "File::Copy", "Parallel::ForkManager", "Bio::SeqIO", "Bio::Seq", "Text::Levenshtein::XS", "Text::WagnerFischer", "List::Util", "Statistics::Basic", "Excel::Writer::XLSX", "Math::Round");
foreach my $module (@modules)
 {
   if (try_load($module)) {
    print "$module is already installed\n";
   } else {
    print "$module is not installed\n";
    CPAN::Shell->install($module);
    #system("ppm install ");
  }
}

sub try_load {
  my $mod = shift;    
  eval("use $mod");

  if ($@) {
        #print "\$@ = $@\n";
        return(0);
  } else {
        return(1);
  }
}


##end##
#
######Author#####
##Jianshi Frank Jin
#
######Version#####
##V1.002
##2022.10.20

