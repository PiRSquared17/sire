#!/bin/env perl

use File::Compare;


## Script to check two directory trees to ensure that they 
## are identical - and to copy files from the first tree
## that are different or missing in the second tree
##
## (C) Christopher Woods, 2014

sub compareDirs
{
      my ($dir1,$dir0) = @_;

      #scan the contents of dir1 and ensure that they
      #all appear in dir0
      local(*DIR1);
      
      opendir DIR1,"$dir1" or die "Cannot open '$dir1': $!";
      
      @files = readdir(DIR1);
      
      foreach $file (@files)
      {
            if ($file =~ m/^\.|~$|\.o$/o)
            {
                  next;
            }

            $fullfile1 = "$dir1/$file";
            $fullfile0 = "$dir0/$file";
            
            if ( ! -e $fullfile0 )
            {
                  print "\n$fullfile1 exists but $fullfile0 doesn't!\n";
                  $cmd = "cp -a $fullfile1 $fullfile0";
                  print "$cmd\n";
                  system($cmd);
            }
            elsif( -d $fullfile0 )
            {
                  #print "\nDirectory $fullfile0...\n";
            }
            elsif( compare($fullfile0,$fullfile1) != 0 )
            {
                  print "\n$fullfile1 exists but is different to $fullfile0!\n";
                  $cmd = "cp $fullfile1 $fullfile0";
                  print "$cmd\n";
                  system($cmd);
            }
      }
}


$dir0 = $ARGV[0];
$dir1 = $ARGV[1];

compareDirs($dir0,$dir1);
